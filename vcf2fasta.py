
from io import StringIO, BytesIO
from tempfile import NamedTemporaryFile, TemporaryDirectory
from subprocess import check_call, Popen, PIPE, check_output
import os
import gzip
from collections import OrderedDict
from itertools import zip_longest
import unittest
from functools import partial
from multiprocessing import Pool

from Bio import SeqIO
from Bio.Data.IUPACData import ambiguous_dna_values


# vcf example alleles modified from:
# http://www.ensembl.org/info/docs/tools/vep/vep_formats.html
#      0 1 2 3 4
# Ref: a t C g a // C is the reference base
# 1  : a t G g a // C base is a G in individual 1
# 2  : a t - g a // C base is deleted w.r.t. the reference in individual 2
# 3  : a t CAgTT // A base is inserted w.r.t. the reference sequence in individual 3
# VCF equivalent
# indi1: 20   3   .   C   G   .   PASS   .
# indi2: 21   2   .   TC   T   .   PASS   .
# indi3: 22   3   .   C   CA   .   PASS   .
# indi3: 22   4   .   g   gT   .   PASS   .
# indi3: 22   5   .   a   T   .   PASS   .

# python 0 1 2 3 4
# VCF    1 2 3 4 5 6 7 8
# Ref:   a c g t a c g t
# 1  :   A - g A C c g -
# 2  :   a c g t A - g a
# 1  :   - -   | |   - -
# 2  :           - -   |


VCF1 = b'''#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2\tind3
20\t3\t.\tC\tG\t.\tPASS\t.\tGT\t1/1\t0/0\t0/.
21\t2\t.\ttC\tt\t.\tPASS\t.\tGT\t0/0\t1/1\t0/0
22\t3\t.\tc\tcA\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1
22\t4\t.\tg\tgT\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1
22\t5\t.\ta\tT\t.\tPASS\t.\tGT\t0/0\t1/1\t.'''

OVERLAPING_SNPS_VCF = b'''#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1
2\t2\t.\tTC\tT\t.\tPASS\t.\tGT\t1/1
2\t4\t.\tC\tG\t.\tPASS\t.\tGT\t1/1
2\t4\t.\tC\tT\t.\tPASS\t.\tGT\t0/0'''

OVERLAPING_SNPS2_VCF = b'''#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti1\ti2
2\t1\t.\tAC\tA\t.\tPASS\t.\tGT\t1/1\t0/0
2\t4\t.\tT\tA\t.\tPASS\t.\tGT\t1/1\t0/0
2\t5\t.\tA\tC\t.\tPASS\t.\tGT\t1/1\t0/0
2\t5\t.\tAC\tA\t.\tPASS\t.\tGT\t0/0\t1/1
2\t7\t.\tGT\tG\t.\tPASS\t.\tGT\t1/1\t0/0
2\t8\t.\tT\tA\t.\tPASS\t.\tGT\t0/0\t1/1
'''

VCF_DP = b'''#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti1\ti2\ti3
20\t3\t.\tC\tG\t.\tPASS\t.\tGT:DP\t1/1:5\t0/0:9\t0/0:1'''


FASTA1 = '''>2\natCga\n>20\natCga\n>21\natCga\n>22\natCga\n'''

FASTA2 = '''>2\nacgtACGT\n'''


# TODO check that the reference seq matches the ref allele
# TODO accept coverages per base per sample to mask sequences


def region_to_str(region):
    if len(region) == 3:
        region = '%s:%d-%d' % (region[0], int(region[1]) + 1, region[2])
    elif len(region) == 2:
        region = '%s:%d' % (region[0], int(region[1]) + 1)
    elif len(region) == 1:
        region = '%s' % (region[0])
    return region


class VCF():
    def __init__(self, vcf_fpath):
        self.vcf_fpath = vcf_fpath
        self.samples = self._read_samples()

    def _read_samples(self):
        for line in gzip.open(self.vcf_fpath):
            if not line.startswith(b'#CHROM'):
                continue
            return line.split()[9:]

    def get_snps_from_region(self, region):
        '''start and end 0-based. end is not included.'''

        region = region_to_str(region)

        tabix = ['tabix', self.vcf_fpath, region]

        stdout_fhand = NamedTemporaryFile(suffix='.vcf')
        process = Popen(tabix, stdout=stdout_fhand)
        process.wait()
        stdout_fhand.seek(0)
        return stdout_fhand, self._parse_vcf(stdout_fhand)

    def _parse_vcf(self, fhand):
        for line in fhand:
            items = line.strip().split(b'\t')
            chrom = items[0]
            start = int(items[1]) - 1
            ref = items[3]
            alt = items[4].split(b',')
            alleles = [ref] + alt

            length = len(ref)
            stop = start + length

            fmt = items[8]
            calls = items[9:]
            yield {'chrom': chrom,
                   'start': start,
                   'stop': stop,
                   'alleles': alleles,
                   'calls': self._parse_calls(fmt, calls)}

    def _parse_calls(self, fmt, calls):
        fmt = fmt.split(b':')
        parsed_calls = []
        for call in calls:
            call = dict((zip(fmt, call.split(b':'))))
            gt = call[b'GT']
            if b'/' in gt:
                gt = [_allele_to_int(allele) for allele in gt.split(b'/')]
            else:
                gt = None
            call[b'GT'] = gt
            parsed_calls.append(call)
        calls = OrderedDict(zip(self.samples, parsed_calls))
        return calls


def _allele_to_int(allele):
    if allele == b'.':
        return None
    else:
        return int(allele)


def parse_bed(fhand):
    for line in fhand:
        items = line.split()
        chrom = items[0]
        start = int(items[1])
        stop = int(items[2])
        yield chrom, start, stop


class VCF2Seq():
    def __init__(self, ref_seqs, vcf, regions):
        self.ref_seqs = ref_seqs
        self.vcf = vcf
        self.regions = regions

    def seqs_per_sample():
        for region in regions:
            snps = self.vcf.get_snps_from_region(region)
            ref_seq = ref_seqs.get_region(region)
            generate_seqs_for_samples(ref_seq, snps, region)


class PeekableIterator(object):
    def __init__(self, iterable):
        self._stream = iterable
        self._buffer = []

    def __iter__(self):
        return self

    def __next__(self):
        if self._buffer:
            item = self._buffer.pop(0)
        else:
            item = self._stream.__next__()
        return item

    def peek(self):
        try:
            item = self._stream.__next__()
        except StopIteration:
            raise
        self._buffer.append(item)
        return item


def _get_overlaping_snps(init_snp, snps):
    pos = init_snp['stop']
    overlaping_snps = [init_snp]
    while True:
        try:
            snp = snps.peek()
        except StopIteration:
            break
        if snp['start'] < pos:
            overlaping_snps.append(snp)
        else:
            break

    for _ in range(len(overlaping_snps) - 1):
        try:
            snp = snps.__next__()
            # print('purging', snp['start'], snp['stop'])
        except StopIteration:
            msg = 'Wrong turn, we should never be here'
            raise RuntimeError(msg)
    return overlaping_snps


def _split_segments_btw_snps(ref_seqs, snps, region):

    snps = PeekableIterator(snps)

    sample_ref_seqs = ref_seqs.get_masked_region(region)

    offset = int(region[1]) if len(region) > 1 else 0
    debug = False
    if debug:
        print('offset', offset)

    pos = offset
    n_samples = ref_seqs.n_samples

    if debug:
        print('*' * 20)
        print('region', region)
    for snp in snps:
        if debug:
            print('pos', pos, 'offset', offset)
            print('snp', snp)
            print('snp_start stop', snp['start'], snp['stop']   )
            print('snp->', pos, snp['start'], snp['stop'])
        if snp['start'] < pos:
            if debug:
                print('overlaping indel at beging')
            # we're at the beginig of a segment, but there's an
            # indel that started before and overlaps with it. We
            # don't yield this snp because is not completely inside
            # the segment
            pos = snp['stop']
        else:
            snps_to_yield = _get_overlaping_snps(snp, snps)
            if snp['start'] == pos:
                if debug:
                    print('segment just after a snp', snp['start'])
                segment_smpl_seqs = [b''] * n_samples
                segment = 0, 0
            else:
                if debug:
                    print('std segment and snp',
                          pos - offset, snp['start'] - offset)
                sstart = pos - offset
                sstop = snp['start'] - offset
                segment = sstart + offset, sstop + offset
                segment_smpl_seqs = [srfs[sstart: sstop] for srfs in sample_ref_seqs]
                #segment_seq = ref_seq[sstart: sstop]
            pos = max([snp_['stop'] for snp_ in snps_to_yield])
            if len(region) < 3 or pos <= region[2]:
                # Are the SNPs inside the region to return?
                if debug:
                    print('yield', segment_smpl_seqs,
                          [(snp_['start'], snp_['stop'] - 1) for snp_ in snps_to_yield])
                yield  segment, segment_smpl_seqs, snps_to_yield

    reamining_smpl_seqs = [srfs[pos - offset:] for srfs in sample_ref_seqs]
    #remaining_seq = ref_seq[pos - offset:]
    if reamining_smpl_seqs[0]:
        reg = pos, region[2] if len(region) > 2 else None
        yield reg, reamining_smpl_seqs, []


IUPAC = {tuple(sorted(nucls)): iupac.encode('utf8') for iupac, nucls in ambiguous_dna_values.items()}


def to_str(bytes_):
    if isinstance(bytes_, bytes):
        return bytes_.decode('utf-8')
    else:
        return bytes_


def to_bytes(str_):
    if not isinstance(str_, bytes):
        return str_.encode('utf-8')
    else:
        return str_


def create_seq_gt(gts, alleles):

    if len(gts) > 1:
        # for the ovelaping snps we just put Ns
        length = max(len(allele) for snp_alleles in alleles for allele in snp_alleles)
        return b'N' * length
    gt = gts[0]
    alleles = alleles[0]

    snp_len = max(len(allele) for allele in alleles)

    if gt is None or None in gt:
        return b'N' * snp_len

    if len(set(gt)) == 1:
        seq = alleles[gt[0]]
    else:
        # hets
        seq = b''
        for nucls_in_pos in zip_longest(*[iter(to_str(alleles[gallele])) for gallele in gt], fillvalue=None):
            if None in nucls_in_pos:
                nucl = b'N'
            else:
                nucls_in_pos = tuple(sorted(set([nuc_ for nuc_ in nucls_in_pos])))
                nucl = IUPAC[nucls_in_pos]
            seq += nucl

    seq = seq.ljust(snp_len, b'-')

    return seq            


def _sum_strs(seq1, seq2):
    try:
        return seq1 + seq2
    except TypeError:
        if isinstance(seq1, bytes):
            return seq1 + seq2.encode("utf-8", "strict")
        else:
            return seq1.encode("utf-8", "strict") + seq2


def _get_gts_and_alleles_for_sample(snps, sample, min_gt_dp):
    gts = []
    alleles = []
    for snp in snps:
        if snp is None:
            gt = None
        else:
            call = snp['calls'][sample]
            if min_gt_dp is not None and int(call[b'DP']) < min_gt_dp:
                gt = None
            else:
                gt = call[b'GT']
        gts.append(gt)
        alleles.append(snp['alleles'])
    return gts, alleles


def generate_seqs_for_samples(region, ref_seqs, vcf, min_gt_dp=None):

    samples = vcf.samples

    temp_file, snps = vcf.get_snps_from_region(region)
    sample_seqs = [b''] * len(samples)
    for sub_region, smpl_segment_seqs, next_snps in _split_segments_btw_snps(ref_seqs,
                                                                   snps,
                                                                   region):
        for isample, sample in enumerate(samples):
            sample_seqs[isample] = _sum_strs(sample_seqs[isample],
                                             smpl_segment_seqs[isample])
            if next_snps:
                gts, alleles = _get_gts_and_alleles_for_sample(next_snps,
                                                               sample,
                                                               min_gt_dp)
                gt = create_seq_gt(gts, alleles)
                sample_seqs[isample] = _sum_strs(sample_seqs[isample], gt)
    return region, list(zip(samples, sample_seqs))


def write_regions_in_fasta(seq_regions, out_dir):
    for region, seqs in seq_regions:
        region = region_to_str(region)
        fpath = os.path.join(out_dir, region + '.fasta')
        fhand = open(fpath, 'wb')
        for indi, indi_seq in seqs:
            fhand.write(b'>')
            fhand.write(to_bytes(region))
            fhand.write(b'\n')
            fhand.write(indi_seq)
            fhand.write(b'\n')
        fhand.close()


def vcf2fasta(vcf_fpath, fasta_fpath, bed_fhand, out_dir, coverages=None,
              min_gt_dp=None, n_threads=None):
    vcf = VCF(vcf_fpath)
    ref_seqs = SeqsSam(fasta_fpath, coverages=coverages,
                       n_samples=len(vcf.samples))
    regions = parse_bed(bed_fhand)
    gen_seqs_for_region = partial(generate_seqs_for_samples, ref_seqs=ref_seqs,
                                  vcf=vcf, min_gt_dp=min_gt_dp)
    if n_threads is None:
        region_seqs = map(gen_seqs_for_region, regions)
    else:
        with Pool(processes=n_threads) as pool:
            region_seqs = pool.imap_unordered(gen_seqs_for_region, regions)
    write_regions_in_fasta(region_seqs, out_dir)


class _SamtoolsSampleCoverages():
    def __init__(self, fhand, sep):
        self._lines = PeekableIterator(fhand)
        self._peeked_items = None
        self.sep = sep
        self.n_samples = self._get_n_samples()

    def _get_n_samples(self):
        items = self._peek_items()
        return len(items[2])

    def _go_to_pos(self, chrom, pos):

        while True:
            peeked_items = self._peek_items()
            relative_pos = self._relative_pos1_pos2(peeked_items[0],
                                                    peeked_items[1],
                                                    chrom, pos)
            if relative_pos == 1:
                self._next_items()
            else:
                break

    def _relative_pos1_pos2(self, chrom1, pos1, chrom2, pos2):
        # None is at the end
        if chrom1 is None:
            return -1
        if chrom2 is None:
            return 1
        if chrom1 > chrom2:
            return -1
        if chrom1 < chrom2:
            return 1
        if pos1 > pos2:
            return -1
        if pos1 < pos2:
            return 1
        return 0

    def _parse_cov_line(self, line):
        items = line.rstrip().split(self.sep)
        chrom = items[0]
        pos = int(items[1]) - 1
        covs = [int(item) for item in items[2:]]
        return chrom, pos, covs

    def _empty_items(self, chrom, pos):
        return chrom, pos, [0] * self.n_samples

    def _peek_items(self):
        if self._peeked_items:
            return self._peeked_items

        try:
            line = self._lines.peek()
            items = self._parse_cov_line(line)
        except StopIteration:
            items = self._empty_items(None, None)
        self._peeked_items = items
        return items

    def _next_items(self):
        self._peeked_items = None
        try:
            line = self._lines.__next__()
            items = self._parse_cov_line(line)
        except StopIteration:
            items = self._empty_items(None, None)
        return items

    def get_data(self, chrom, pos):
        chrom = to_bytes(chrom)
        peeked_items = self._peek_items()
        relative_pos = self._relative_pos1_pos2(chrom, pos, peeked_items[0],
                                                peeked_items[1])
        if relative_pos == -1:
            self._go_to_pos(chrom, pos)

        peeked_items = self._peek_items()
        relative_pos = self._relative_pos1_pos2(chrom, pos, peeked_items[0],
                                                peeked_items[1])
        if relative_pos == 0:
            return self._next_items()[2]
        elif relative_pos == 1:
            return self._empty_items(chrom, pos)[2]
        else:
            msg = 'We shold not be here, fixme'
            raise RuntimeError(msg)            


class GenomicCoverages():
    def __init__(self, fhand, min_cov, sep=b'\t'):
        self._genome_covs = _SamtoolsSampleCoverages(fhand, sep)
        self.min_cov = min_cov

    def get_regions(self, region):
        min_cov = self.min_cov
        chrom = region[0]
        start = region[1]
        stop = region[2]

        cov_state = None
        region_start = None
        region_stop = None
        for pos in range(region[1], region[2]):
            covs_in_pos = self._genome_covs.get_data(chrom, pos)
            covs_in_pos = [True if cov >= min_cov else False for cov in covs_in_pos]
            if cov_state is None:
                cov_state = covs_in_pos
                region_start = pos
                continue
            if cov_state == covs_in_pos:
                continue
            else:
                yield (region_start, pos), cov_state
                cov_state = covs_in_pos
                region_start = pos
        else:
            if region_start != pos + 1:
                yield (region_start, pos + 1), cov_state


class Seqs():
    def __init__(self, coverages=None, n_samples=None):
        self.coverages = coverages
        self.n_samples = n_samples

    def get_masked_region(self, region):
        '''It returns the masked sequence for every sample.

        The masking is done according to the BAM coverage.'''
        seq = self.get_region(region)

        if len(region) < 3:
            chrom = region[0]
            start = region[1] if len(region) > 1 else 0
            end = region[2] if len(region) > 2 else len(seq) + start
            region = [chrom, start, end]
        else:
            chrom, start, end = region

        n_samples = self.n_samples

        if self.coverages is None:
            subregions = [((region[1], region[2]), [True] * n_samples)]
        else:
            subregions = self.coverages.get_regions(region)

        masked_seqs = [b''] * n_samples
        for subregion, are_covered in subregions:
            assert len(are_covered) == n_samples

            subseq = seq[subregion[0] - start: subregion[1] - start]
            subseq_len = len(subseq)
            sample_subseqs = [b''] * len(are_covered)
            for sample_idx, sample_is_covered in enumerate(are_covered):
                if sample_is_covered:
                    sample_subseq = subseq
                else:
                    sample_subseq = b'N' * subseq_len
                masked_seqs[sample_idx] += sample_subseq
        return masked_seqs


class SeqsFasta(Seqs):
    def __init__(self, fhand, coverages=None, n_samples=None):
        self.seqs = SeqIO.to_dict(SeqIO.parse(fhand, 'fasta'))
        super().__init__(coverages=coverages, n_samples=n_samples)

    def get_region(self, region):
        chrom = region[0]
        start = int(region[1]) if len(region) > 1 else None
        stop = int(region[2]) if len(region) > 2 else None
        slc = slice(start, stop)

        seq = self.seqs[chrom]
        return str(seq[slc].seq).encode('utf-8')


def parse_fasta(fhand):
    'it iterates over a fasta fhand and yields the content of each fasta item'
    seq         = []
    name        = None
    description = None
    for line in fhand:
        line = line.strip()
        if not line or line[0] == '#':
            continue
        if line.startswith(b'>'):
            if seq:
                seq = b''.join(seq)
                yield name, description, seq
                seq         = []
                name        = None
                description = None

            items = line.split(b' ', 1)
            name  = items[0][1:]
            try:
                description = items[1]
            except IndexError:
                description = None
        else:
            seq.append(line)
    else:
        seq = b''.join(seq)
        yield name, description, seq


def first(iterator):
    for item in iterator:
        return item
    raise ValueError('iterator was empty')


class SeqsSam(Seqs):
    def __init__(self, fpath, coverages=None, n_samples=None):
        self.fpath = fpath
        super().__init__(coverages=coverages, n_samples=n_samples)

    def get_region(self, region):
        sam = ['samtools', 'faidx', self.fpath, region_to_str(region)]
        fhand = NamedTemporaryFile(suffix='.fasta')
        check_call(sam, stdout=fhand)
        fhand.flush()
        fhand.seek(0)
        return first(parse_fasta(fhand))[2]


class TestVcf():
    def __init__(self, vcf_string):
        self.vcf_string = vcf_string
        

    def __enter__(self):
        self._create_vcf()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        fpaths = [self._vcf_fhand.name,
                  self.vcf_fpath,
                  self.vcf_fpath + '.tbi']
        for fpath in fpaths:
            if os.path.exists(fpath):
                os.remove(fpath)

    def _create_vcf(self):
        self._vcf_fhand = NamedTemporaryFile(suffix='.vcf', delete=False)
        self._vcf_fhand.write(self.vcf_string)
        self._vcf_fhand.flush()
        bgzip = ['bgzip', self._vcf_fhand.name]
        check_call(bgzip)
        self.vcf_fpath = self._vcf_fhand.name + '.gz'
        tabix = ['tabix', '-p', 'vcf', self.vcf_fpath]
        check_call(tabix)


class TestFastaFaidx():
    def __init__(self, fasta_string):
        self.fasta_string = fasta_string
        

    def __enter__(self):
        self._create_fasta()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        fpaths = [self.fasta_fhand.name,
                  self.fasta_fhand.name + '.fai']
        for fpath in fpaths:
            if os.path.exists(fpath):
                os.remove(fpath)

    def _create_fasta(self):
        self.fasta_fhand = NamedTemporaryFile(suffix='.fasta', delete=False,
                                              mode='wt')
        self.fasta_fhand.write(self.fasta_string)
        self.fasta_fhand.flush()
        faidx = ['samtools', 'faidx', self.fasta_fhand.name]
        check_call(faidx)


class Test(unittest.TestCase):

    def test_vcf(self):
        with TestVcf(VCF1) as test_vcf:
            vcf = VCF(test_vcf.vcf_fpath)
            temp_file, snps = vcf.get_snps_from_region((20, 2, 3))
            snps = list(snps)
            assert len(snps) == 1
            snp = snps[0]
            assert snp['start'] == 2
            assert snp['stop'] == 3
            assert snp['calls'][b'ind1'][b'GT'] == [1, 1]
            assert snp['alleles'] == [b'C', b'G']
            temp_file, snps = vcf.get_snps_from_region((20, 1, 2))
            assert not list(snps)
            temp_file, snps = vcf.get_snps_from_region((20, 3, 4))
            assert not list(snps)

        with TestVcf(VCF_DP) as test_vcf:
            vcf = VCF(test_vcf.vcf_fpath)
            temp_file, snps = vcf.get_snps_from_region((20, 2, 3))
            snps = list(snps)
            assert len(snps) == 1
            assert first(first(snps)['calls'].values())[b'DP'] == b'5'

    def test_bed(self):
        bed = 'chrom1\t0\t10\nchrom2\t10\t21'
        bed_fhand = StringIO(bed)
        regions = list(parse_bed(bed_fhand))
        assert regions == [('chrom1', 0, 10), ('chrom2', 10, 21)]

    def test_seqs(self):
        ref_fhand = StringIO(FASTA1)
        seqs = SeqsFasta(ref_fhand)
        assert seqs.get_region(['20']) == b'atCga'
        assert seqs.get_region(['21', 1]) == b'atCga'[1:]
        assert seqs.get_region(['22', 1, 3]) == b'atCga'[1:3]

        with TestFastaFaidx(FASTA1) as fasta:
            seqs = SeqsSam(fasta.fasta_fhand.name)
            assert seqs.get_region(['20']) == b'atCga'
            assert seqs.get_region(['21', 1]) == b'atCga'[1:]
            assert seqs.get_region(['22', 1, 3]) == b'atCga'[1:3]

    def test_masked_seq(self):
        # 2  atCga
        # 20 atCga
        # 21 atCga
        cov = b'''20 2 9 9 9
20 3 9 9 9
20 4 6 6 7
'''

        fhand = BytesIO(cov)
        covs = GenomicCoverages(fhand, min_cov=7, sep=b' ')
        ref_fhand = StringIO(FASTA1)
        seqs = SeqsFasta(ref_fhand, coverages=covs, n_samples=3)
        expected = [b'NtCNN', b'NtCNN', b'NtCgN']
        assert list(seqs.get_masked_region(['20'])) == expected

        fhand = BytesIO(cov)
        covs = GenomicCoverages(fhand, min_cov=6, sep=b' ')
        ref_fhand = StringIO(FASTA1)
        seqs = SeqsFasta(ref_fhand, coverages=covs, n_samples=3)
        expected = [b'NNNNN', b'NNNNN', b'NNNNN']
        assert list(seqs.get_masked_region(['2'])) == expected
        expected = [b'tCgN', b'tCgN', b'tCgN']
        assert list(seqs.get_masked_region(['20', 1])) == expected
        expected = [b'NNNNN', b'NNNNN', b'NNNNN']
        assert list(seqs.get_masked_region(['2'])) == expected

        ref_fhand = StringIO(FASTA1)
        seqs = SeqsFasta(ref_fhand, n_samples=3)
        expected = [b'atCga', b'atCga', b'atCga']
        assert list(seqs.get_masked_region(['20'])) == expected

    def test_seq_for_sample(self):
        ref_fhand = StringIO(FASTA1)

        with TestVcf(VCF1) as test_vcf:
            vcf = VCF(test_vcf.vcf_fpath)
            refs = SeqsFasta(ref_fhand, n_samples=len(vcf.samples))
            regions=[['20'], ('21', 1, 5), ['22']]
            generate_seqs_for_region = partial(generate_seqs_for_samples,
                                               ref_seqs=refs, vcf=vcf)
            res = map(generate_seqs_for_region, regions)
            res = list(res)
            assert res == [(['20'], [(b'ind1', b'atGga'),
                                     (b'ind2', b'atCga'),
                                     (b'ind3', b'atNga')]),
                           (('21', 1, 5), [(b'ind1', b'tCga'),
                                           (b'ind2', b't-ga'),
                                           (b'ind3', b'tCga')]),
                            (['22'], [(b'ind1', b'atc-g-a'),
                                      (b'ind2', b'atc-g-T'),
                                      (b'ind3', b'atcAgTN')])]

        with TestVcf(OVERLAPING_SNPS_VCF) as test_vcf:
            vcf = VCF(test_vcf.vcf_fpath)
            generate_seqs_for_region = partial(generate_seqs_for_samples,
                                               ref_seqs=refs, vcf=vcf)
            res = map(generate_seqs_for_region, [['2']])
            try:
                list(res)
                self.fail('RuntimeError expected')
            except:
                pass

        # With DP
        #20 3   .   C   G   .   PASS    .   GT:DP   1/1:5   0/0:9   0/0:1
        with TestVcf(VCF_DP) as test_vcf:
            vcf = VCF(test_vcf.vcf_fpath)
            regions=[['20']]
            generate_seqs_for_region = partial(generate_seqs_for_samples,
                                               ref_seqs=refs, vcf=vcf,
                                               min_gt_dp=3)
            res = map(generate_seqs_for_region, regions)
            res = list(res)
            assert res == [(['20'], [(b'i1', b'atGga'),
                                     (b'i2', b'atCga'),
                                     (b'i3', b'atNga')])]

            generate_seqs_for_region = partial(generate_seqs_for_samples,
                                               ref_seqs=refs, vcf=vcf,
                                               min_gt_dp=10)
            res = map(generate_seqs_for_region, regions)
            res = list(res)
            assert res == [(['20'], [(b'i1', b'atNga'),
                                     (b'i2', b'atNga'),
                                     (b'i3', b'atNga')])]


            generate_seqs_for_region = partial(generate_seqs_for_samples,
                                               ref_seqs=refs, vcf=vcf,
                                               min_gt_dp=7)
            res = map(generate_seqs_for_region, regions)
            res = list(res)
            assert res == [(['20'], [(b'i1', b'atNga'),
                                     (b'i2', b'atCga'),
                                     (b'i3', b'atNga')])]

            # 20 atCga
            cov = b'''20 2 9 9 9
20 3 9 9 9
20 4 6 6 6
20 5 6 6 6
'''

            fhand = BytesIO(cov)
            covs = GenomicCoverages(fhand, min_cov=1, sep=b' ')
            ref_fhand = StringIO(FASTA1)
            refs = SeqsFasta(ref_fhand, coverages=covs,
                             n_samples=len(vcf.samples))
            generate_seqs_for_region = partial(generate_seqs_for_samples,
                                               ref_seqs=refs, vcf=vcf,
                                               min_gt_dp=7)
            res = map(generate_seqs_for_region, regions)
            res = list(res)
            assert res == [(['20'], [(b'i1', b'NtNga'),
                                     (b'i2', b'NtCga'),
                                     (b'i3', b'NtNga')])]


            fhand = BytesIO(cov)
            covs = GenomicCoverages(fhand, min_cov=9, sep=b' ')
            ref_fhand = StringIO(FASTA1)
            refs = SeqsFasta(ref_fhand, coverages=covs,
                             n_samples=len(vcf.samples))
            generate_seqs_for_region = partial(generate_seqs_for_samples,
                                               ref_seqs=refs, vcf=vcf,
                                               min_gt_dp=9)
            res = map(generate_seqs_for_region, regions)
            res = list(res)
            assert res == [(['20'], [(b'i1', b'NtNNN'),
                                     (b'i2', b'NtCNN'),
                                     (b'i3', b'NtNNN')])]

            fhand = BytesIO(cov)
            covs = GenomicCoverages(fhand, min_cov=10, sep=b' ')
            ref_fhand = StringIO(FASTA1)
            refs = SeqsFasta(ref_fhand, coverages=covs,
                             n_samples=len(vcf.samples))
            generate_seqs_for_region = partial(generate_seqs_for_samples,
                                               ref_seqs=refs, vcf=vcf,
                                               min_gt_dp=7)
            res = map(generate_seqs_for_region, regions)
            res = list(res)
            assert res == [(['20'], [(b'i1', b'NNNNN'),
                                     (b'i2', b'NNCNN'),
                                     (b'i3', b'NNNNN')])]


    def _get_segments_start_end(self, segments):
        res = []
        for reg, seq, snps in segments:
            start_stops = []
            if snps:
                start_stops = [(snp['start'], snp['stop'] - 1) for snp in snps]
            res.append((reg, seq, start_stops))
        return res

    def test_split_segments(self):
        #          01234567
        # FASTA2 2 acgtACGT

        ref_fhand = StringIO(FASTA2)
        seqs = SeqsFasta(ref_fhand, n_samples=1)
        with TestVcf(OVERLAPING_SNPS2_VCF) as test_vcf:
            vcf = VCF(test_vcf.vcf_fpath)
            regions = [('2',)]
            expected = [[((0, 0), [b''], [(0, 1)]),
                         ((2, 3), [b'g'], [(3, 3)]),
                         ((0, 0), [b''], [(4, 4), (4, 5)]),
                         ((0, 0), [b''], [(6, 7), (7, 7)])]]

            regions.append(('2', 1))
            expected.append([((2, 3), [b'g'], [(3, 3)]),
                             ((0, 0), [b''], [(4, 4), (4, 5)]),
                             ((0, 0), [b''], [(6, 7), (7, 7)])])

            regions.append(('2', 1, 4))
            expected.append([((2, 3), [b'g'], [(3, 3)])])

            regions.append(('2', 1, 5))
            expected.append([((2, 3), [b'g'], [(3, 3)])])

            for region, exp in zip(regions, expected):
                temp_file, snps = vcf.get_snps_from_region(region)
                segments = list(_split_segments_btw_snps(seqs, snps, region))
                assert self._get_segments_start_end(segments) == exp

    def test_seq_gt(self):
        
        assert create_seq_gt([[0, 0]], [[b'A', b'T']]) == b'A'
        assert create_seq_gt([[1, 1]], [[b'A', b'T']]) == b'T'
        assert create_seq_gt([[0, 1]], [[b'A', b'T']]) == b'W'
        assert create_seq_gt([None], [[b'A', b'T']]) == b'N'
        assert create_seq_gt([[0, None]], [[b'A', b'T']]) == b'N'
        assert create_seq_gt([None], [[b'A', b'ATT']]) == b'NNN'
        assert create_seq_gt([[0, 1]], [[b'A', b'ATT']]) == b'ANN'
        assert create_seq_gt([[0, 0]], [[b'A', b'ATT']]) == b'A--'
        assert create_seq_gt([[0, 1]], [[b'ATT', b'A']]) == b'ANN'
        assert create_seq_gt([[0, 0], [0, 0]],
                              [[b'ATT', b'A'], [b'A', b'T']]) == b'NNN'

    def test_write_regions(self):
        res = [(['20'], [(b'ind1', b'atGga'),
                         (b'ind2', b'atCga'),
                         (b'ind3', b'atNga')]),
               (('21', 1, 5), [(b'ind1', b'tCga'),
                               (b'ind2', b'tga'),
                               (b'ind3', b'tCga')]),
               (['22'], [(b'ind1', b'atcga'),
                         (b'ind2', b'atcgT'),
                         (b'ind3', b'atcAgTN')])]

        with TemporaryDirectory() as out_dir:
            write_regions_in_fasta(res, out_dir)
            fpath = os.path.join(out_dir, first(sorted(os.listdir(out_dir))))
            fhand = open(fpath)
            assert fhand.read() == '>20\natGga\n>20\natCga\n>20\natNga\n'
            fhand.close()

    def test_main(self):
        with TestVcf(VCF1) as test_vcf:
            with TestFastaFaidx(FASTA1) as fasta:
                with TemporaryDirectory() as out_dir:
                    bed = '2\t0\t10\n20\t10\t21'
                    bed_fhand = StringIO(bed)
                    vcf2fasta(test_vcf.vcf_fpath, fasta.fasta_fhand.name,
                              bed_fhand, out_dir)
        return
        with TestVcf(VCF1) as test_vcf:
            with TestFastaFaidx(FASTA1) as fasta:
                with TemporaryDirectory() as out_dir:
                    bed = '2\t0\t10\n20\t10\t21'
                    bed_fhand = StringIO(bed)
                    vcf2fasta(test_vcf.vcf_fpath, fasta.fasta_fhand.name,
                              bed_fhand, out_dir, n_threads=2)

    def test_bam_cov(self):

        cov = b'''20 2 9 9 9
20 3 9 9 9
20 4 6 6 7
'''
        fhand = BytesIO(cov)
        covs = GenomicCoverages(fhand, min_cov=7, sep=b' ')
        expected = [((1, 3), [True, True, True]),
                    ((3, 4), [False, False, True]),
                    ((4, 5), [False, False, False])]
        assert list(covs.get_regions((b'20', 1, 5))) == expected

#chr1 49 1 1 1
#chr1 50 1 1 1
#chr1 51 6 6 6
#chr1 52 6 6 0
#chr1 54 6 6 6
#chr1 55 6 6 7
#chr2 56 6 6 6

        cov = b'''chr1 50 1 1 1
chr1 51 1 1 1
chr1 52 6 6 6
chr1 53 6 6 0
chr1 55 6 6 6
chr1 56 6 6 7
chr1 57 6 6 6
chr2 58 6 6 6
'''
        fhand = BytesIO(cov)
        # region not covered
        covs = _SamtoolsSampleCoverages(fhand, sep=b' ')
        assert covs.get_data(b'chr1', 40) == [0, 0, 0]
        assert covs.get_data(b'chr1', 50) == [1, 1, 1]
        assert covs.get_data(b'chr1', 51) == [6, 6, 6]
        assert covs.get_data(b'chr1', 53) == [0, 0, 0]
        assert covs.get_data(b'chr1', 55) == [6, 6, 7]
        assert covs.get_data(b'chr1', 56) == [6, 6, 6]
        assert covs.get_data(b'chr1', 70) == [0, 0, 0]
        assert covs.get_data(b'chr2', 70) == [0, 0, 0]

        fhand = BytesIO(cov)
        covs = GenomicCoverages(fhand, min_cov=1, sep=b' ')
        expected = [((48, 49), [False, False, False]),
                    ((49, 52), [True, True, True])]
        assert list(covs.get_regions((b'chr1', 48, 52))) == expected

        fhand = BytesIO(cov)
        covs = GenomicCoverages(fhand, min_cov=1, sep=b' ')
        expected = [((40, 45), [False, False, False])]
        assert list(covs.get_regions((b'chr1', 40, 45))) == expected

        fhand = BytesIO(cov)
        covs = GenomicCoverages(fhand, min_cov=1, sep=b' ')
        expected = [((48, 49), [False, False, False]),
                    ((49, 52), [True, True, True]),
                    ((52, 53), [True, True, False]),
                    ((53, 54), [False, False, False]),
                    ((54, 57), [True, True, True]),
                    ((57, 60), [False, False, False])]
        assert list(covs.get_regions((b'chr1', 48, 60))) == expected

        fhand = BytesIO(cov)
        covs = GenomicCoverages(fhand, min_cov=1, sep=b' ')
        expected = [((49, 52), [True, True, True]),
                    ((52, 53), [True, True, False]),
                    ((53, 54), [False, False, False]),
                    ((54, 57), [True, True, True])]
        assert list(covs.get_regions((b'chr1', 49, 57))) == expected


def build_tomato_phylome_genomes():
    out_dir = 'acc_build_genomes'
    vcf_fpath = 'phylome.20161002.dp3.sample.vcf.gz'
    ref_genome_fpath = '/home/jope/genomes/tomato/S_lycopersicum_chromosomes.2.50.fa'
    bed_fpath = '/home2/jope/tomato/snv_calling/tomato_cov_analisys/region_cov_less_200and_complex.bed'
    bed_fhand = open(bed_fpath)
    threads = None
    vcf2fasta(vcf_fpath, ref_genome_fpath, bed_fhand, out_dir, min_gt_dp=3,
              n_threads=threads)

if __name__ == '__main__':
    #import sys; sys.argv = ['', 'Test.test_seq_for_sample']
    unittest.main()
    #build_tomato_phylome_genomes()

