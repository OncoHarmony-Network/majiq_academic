Splicegraph Modules
===================

Splicegraph modules (Modules) were first conceptualized and applied to
analyzing RNA-Seq data by Hu et al (1). Modules are derived as follows:
as you trace a gene from 5’ to 3’, a Module begins when there are two or
more possible ways to splice out of a locus, and the Module ends when
the alternative splicing paths converge to a single locus. For example,
a simple cassette exon represented by a Module would start at the 5’
splice site of the upstream constitutive exon, and the Module would end
at the 3’ splice site of the downstream constitutive exon.

We applied the Module logic to identify simple alternative
isoform-forming events (two possible paths in the Module: cassette
exons, tandem cassette exons, binary alternative splice sites, mutually
exclusive exons, intron retentions, alternative first/last exons, etc).
Furthermore, the Module logic is easily extendable to capture more
complex alternative splicing events (three or more possible paths in the
Module).

The first step of the VOILA Categorizer is to identify modules along
every gene. Each module has an ID which is simply: GeneID_<N> where <N>
is the nth module found while scanning the gene in the 5’ to 3’
direction.

The second step of the voila Categorizer is to identify what types of
alternative-isoform-generating events are within each module.

