#############
MAJIQ Builder
#############

.. role:: bash(code)
   :language: bash

The MAJIQ Builder defines a splicegraph and coverage over local splicing
variations per input experiment.
The MAJIQ Builder pipeline can be run on a single machine using the command
:bash:`majiq build`.
This command chains together the following steps:

1. Initialize splicegraph using annotated transcript definitions from GFF3
   (:bash:`majiq-build gff3`).
2. Process input BAM files for spliced junction and retained intron coverage
   (:bash:`majiq-build sj`).
3. Update splicegraph, identifying annotated/novel junctions/retained introns
   with reliable coverage (:bash:`majiq-build update`).
4. Get raw and bootstrapped coverage per LSV, experiment for downstream
   quantification/analysis (:bash:`majiq-build psi-coverage`).


Input/output files
==================


From the perspective of input and output files, the MAJIQ builder takes as
input and produces as output the following:

Input files
-----------

- GFF3 format annotation database defining annotated transcripts. This defines
  the annotated exons, junctions, and retained introns which initialize the
  output splicegraph.
- Input RNA-seq experiments to obtain coverage over LSVs for. This can either
  be in the form of BAM files (with contigs matching the GFF3 file) or "SJ"
  files produced by previous runs of the MAJIQ builder using the same GFF3.

Output files
------------

- Splicegraph database (`splicegraph.zarr`) defining annotated and novel exons,
  junctions, and retained introns on which LSVs are defined.
- SJ files (`{experiment}.sj`) reusable intermediate files with coverage
  information per experiment/BAM.
  Used as input by subsequent runs of the MAJIQ builder (used to speed up
  execution of future builds with the same experiment and GFF3 by parsing input
  BAM file only once).
- PsiCoverage files (`{group}.psicov`) with raw/bootstrapped coverage over LSVs
  for one or more experiments saved as a group.
  Used as input by MAJIQ quantifiers and other downstream analysis.


MAJIQ builder pipeline CLI
==========================

From the perspective of the command-line, MAJIQ specifies the locations of
these inputs and outputs using the required arguments:

- :bash:`gff3`: path to input GFF3 file.
- :bash:`experiments_tsv`: path to TSV file defining where to find input
  experiments and how they will be grouped together for building the
  splicegraph (see below for details).
- :bash:`output_dir`: path for new output directory to which splicegraph, SJ,
  PsiCoverage files will be saved.

The `experiments_tsv` has required column `path` and optional columns
`group` and `strandness`:

- The `path` column indicates the input BAM/SJ file for each input experiment.
- The `group` column defines names of independent "build groups" that MAJIQ
  processes together.
  If omitted, all experiments are treated as a single build group.
- The `strandness` column, when included and specified for an experiment,
  overrides the :bash:`--strandness` flag controlling parsing of
  strand-specific coverage from BAM files.

So an example over five experiments, explicitly specifying strandness for two
experiments (A1, A3) and reusing coverage from two other experiments (B1, C2)
could be:

====== ===== ==========
path   group strandness
====== ===== ==========
A1.bam A     REVERSE
A2.bam A
A3.bam A     FORWARD
B1.sj  B
C2.sj  C
====== ===== ==========

Build groups are used when multiple experiments with evidence for a junction or
retained intron should be required before they are considered reliable (whether
annotated or *de novo*).
The set of reliable junctions and retained introns determines the splicegraph
and which LSVs are saved to output majiq files for quantification.
When analyzing patterns summarizing groups of experiments (replicates of a cell
line/condition or samples from the same tissue type), grouping them together is
often appropriate.
This allows focus on evidence found in multiple experiments.
However, when analyzing variation between individual samples (no replicates,
differences between samples from the same tissue type), grouping samples
independently may be more appropriate.
This allows analysis of changes found in single experiments.
The most important optional parameter governing analysis of these build groups
is :bash:`--min-experiments`, which specifies how many (or what proportion) of
experiments in a build group are required to provide evidence for a reliable
intron or junction.
Note when :bash:`--min-experiments 1` that there is no difference between
grouping experiments together vs independently, as a single experiment from any
build group will then provide sufficient evidence.

We believe our defaults are sensible, but it is worth paying particular
attention to the following parameters:

- :bash:`--min-experiments`: as explained above
- :bash:`--mindenovo`: minimum readrate to pass a novel junction or retained
  intron into the splicegraph
- :bash:`--simplify`: ignore reliable but very low usage junctions or retained
  introns
- :bash:`--quantify-minreads`: minimum readrate for a junction or retained
  intron to pass an LSV for quantification

More detailed explanations of these parameters (and others) can be found by
running :bash:`majiq build --help`.


Finer control with :bash:`majiq-build`
======================================

The MAJIQ Builder pipeline chains together 4 different unique commands
from :bash:`majiq-build`.
There are many cases where you might want to use these commands directly rather
than the pipeline.
These include (but are not limited to):

- Processing GFF3 and SJ files:
  The initial steps for processing annotations and input BAM files can
  generally be shared between analyses.
  Furthermore, each BAM file can be processed in parallel in a
  cluster/distributed environment.
- Creating PsiCoverage files faster:
  Each PsiCoverage file can be processed in parallel in a cluster/distributed
  environment.
  The pipeline saves experiments with each other in a single file per build
  group.
  Beyond parallelizing over these groups, we observe further speed improvements
  on large groups by splitting them into smaller batches.
- Different group definitions for PsiCoverage files:
  We sometimes want to group experiments for quantification differently than the
  build groups used for updating the splicegraph.
- PsiCoverage from experiments that were not part of a previous build:
  The pipeline only creates coverage for experiments contributing to the
  splicegraph at the time of the build.
  New experiments can be quickly compared to previous quantifications by
  producing coverage relative to an old splicegraph.

The last point is of particular interest.
While :bash:`majiq-build psi-coverage` allows creating PsiCoverage for new
experiments relative to old splicegraphs, by itself it does not evaluate what
the analysis would be if the build had included the additional experiment.


Two-pass build
--------------

The MAJIQ builder includes tools for contrasting and combining multiple
splicegraphs:

- :bash:`majiq-build combine` allows combining independent evidence from
  multiple splicegraphs into a single splicegraph.
  This is roughly equivalent to running :bash:`majiq-build update` with the
  experiments from each build as independent build groups.
  Simplification is not exactly the same.
  **NOTE**: how we handle introns needs to be fixed. We can propagate intron
  status to all introns between annotated exons to enable near-equivalence
  (besides slight differences with simplification) as a fix.
- :bash:`majiq-build combine` allows treating novel junctions from some of
  these splicegraphs as known, highlighting junctions that were novel to
  specific experiments.
- :bash:`majiq-build psi-coverage` allows producing coverage for events that
  are unique to only one splicegraph (i.e. if it was structurally the same in
  the first build, ignore it).
  This enables focusing on structurally novel events.
  It can also prevent duplicate work with shared experiments/events which were
  quantified in previous builds that share the same events.

This functionality is of particular interest for our clinical analysis
pipelines for patients with suspected Mendelian disorders, where each
per-patient analysis shares a large group of controls.
The controls can be analyzed one time (first pass), with a second pass analysis
for each patient afterwards.
