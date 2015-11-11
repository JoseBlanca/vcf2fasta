vcf2fasta creates sequences for each sample in a VCF file.

I have written this script to help me with a phylome analysis.

It requires:

    * a genome reference fasta file
    * a VCF file.
    * a BED file with the regions to process

It will substitute the alleles found in the VCF file into the genome reference to create the sample sequences.

It will create a fasta file for each region found in the BED file. In each fasta file the aligned sequence for each sample in the VCF file will be written.

A file with the coverages found in a list of BAM files can be optionally given.
This file can be generated with the command:

    $ samtools depth -f bam_list.txt -q 20 -Q 55 | gzip > coverages.csv.gz

If this coverage file is given all the positions with a coverage lower than a given threshold will be set to N.

Other alternatives to this script are: GATK's [FastaAlternateReferenceMaker](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_fasta_FastaAlternateReferenceMaker.php) and [bcftools consensus](https://samtools.github.io/bcftools/bcftools.html#consensus).

