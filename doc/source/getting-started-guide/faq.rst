.. _faq:

Frequently asked questions
==========================

Why is majiq builder returning 0 or a very small number of LSVs?
    The reason for a very small number of LSVs can be different. The most common one is that annotation DB and the bam files does not match. The bam files can be mapped with a different DB version or the bam files has different names or variation of the contig id. The most common in mammals is that your DB identifies the chromosomes as 1,2,3,4... but the bam identifies them as chr1, chr2, chr3... The way to know if this is your case you can run samtools using the following command, samtools view bamfile | grep -n 200 and see the contig/chromosome ids in the 3rd column. Those names should match with the ones in the annotation DB in the gff3 file.

When running deltapsi/psi/het, why do I get an error like "LSV type does not match"
    MAJIQ analizes the bam files finding de novo junctions and retain introns. That makes dificult its quantification if not all the lsvs are defined in the same way. For that reason we need to build together all the samples that we want to compare, so they start with a common lsv definitions

Does MAJIQ accept different sets of RNA-Seq with different read lengths of strandness types?
    Yes, for the read length you only need to specify in the settings file, the longest readlength in the set of samples. For the strandness type check the majiq quick guide to check how to do it.

What is a half-exon?
    MAJIQ detects de novo junctions that can create new exons. Sometimes the boundaries of those exons are not clear ( not present, to far away, or not well covered...), in those coses we define the exon as half-exon. The missing coordinate is specified as na

In VOILA, some of the LSV type thumbnails don't have a numver, only "unk". What does this mean?
    "unk" stands for unknown. That reflects an lsv where its reference exon is a half-exon.

What difference is there in running PSI on individual samples vs. as a group of replicates for a given condition/cell type?
    MAJIQ PSI and deltapsi can work with one or more replicas per group. The number of replicas makes an effect in the quantification in two main aspects. First the quantification filters out any LSV that does not pass the filter in set of experiments ( The number of experiments is defined by min_emperiments argument, or 50% of the experiments by default). The second aspect where the replicas makes an effect is on the confidence. When MAJIQ quantifies replicas it aggregates the evidence of the replicas in order to improve the confidence on the quantification.




