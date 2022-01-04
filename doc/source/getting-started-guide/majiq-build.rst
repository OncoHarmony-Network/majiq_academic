.. currentmodule:: new_majiq
.. _majiq-build:

##################################
Building splicegraphs and coverage
##################################

.. role:: bash(code)
   :language: bash

The MAJIQ Builder defines a splicegraph over local splicing variations.
The MAJIQ Builder pipeline can be run on a single machine using the command
:bash:`majiq build`.
This command chains together the following steps:

1. Initialize splicegraph using annotated transcript definitions from GFF3
   (:bash:`majiq-build gff3`).
2. Process input BAM files for spliced junction and retained intron coverage
   (:bash:`majiq-build sj`), producing :py:class:`SJExperiment` files.
3. Update splicegraph, identifying annotated/novel junctions/retained introns
   with reliable coverage (:bash:`majiq-build update`),
   producing a :py:class:`SpliceGraph` file.

An additional command is required before quantification.
The command :bash:`majiq psi-coverage` gets raw and bootstrapped coverage per
LSV and experiment for downstream quantification/analysis (:bash:`majiq-build
psi-coverage`) as a :py:class:`PsiCoverage` file.


The MAJIQ Builder pipeline
==========================

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

- :py:class:`SpliceGraph` file (`splicegraph.zarr`) defining
  annotated and novel exons, junctions, and retained introns on which LSVs are
  defined.
- :py:class:`SJExperiment` files (`{experiment}.sj`) reusable
  intermediate files with coverage information per experiment/BAM.
  Used as input by subsequent runs of the MAJIQ builder (used to speed up
  execution of future builds with the same experiment and GFF3 by parsing input
  BAM file only once).


Command-line interface
----------------------

From the perspective of the command-line, the MAJIQ builder pipeline specifies
the locations of these inputs and outputs using the required positional
arguments:

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

More detailed explanations of these parameters (and others) can be found by
running :bash:`majiq build --help`.


.. _majiq-psi-coverage:

The MAJIQ PsiCoverage command
=============================

MAJIQ now requires running :bash:`majiq psi-coverage` separately from the MAJIQ
Builder pipeline.
This makes the builder have the singular purpose of building a splicegraph
model of all genes and leaves preparing coverage for quantification to this new
step.
This also allows for creating PsiCoverage files faster in a cluster/distributed
environment.

This command produces a single PsiCoverage file, requiring a splicegraph and
one or more input SJ files as input:

- :bash:`splicegraph`: path to input :py:class:`SpliceGraph` file
- :bash:`psi_coverage`: path for output :py:class:`PsiCoverage` file
- :bash:`sj [sj ...]`: path for input :py:class:`SJExperiment` files.

We believe our default parameters are sensible, but it is worth paying
particular attention to :bash:`--minreads` and :bash:`--minbins`, which
determine which LSVs pass quantification thresholds. Detailed explanations of
all parameters can be found by running :bash:`majiq psi-coverage --help`.

Coverage from each input experiment (SJ file) is stored independently in the
resulting PsiCoverage file.
Downstream analysis steps permit multiple PsiCoverage files to be grouped
together later.
Since one can group experiments later, creating a PsiCoverage file for each
experiment generally gives the most flexibility for downstream analysis.
However, if the experiments are generally analyzed together as a group,
grouping experiments into smaller number of files allows one (and workflow
managers) to keep track of less commands/output files.
This trades off against distributed parallelism
(since :bash:`majiq psi-coverage` runs on a single machine),
so in these cases it is advisable to split these groups into independent
batches that can be combined later.
For example, in past analyses with GTEx, where some tissues have around 1000
experiments, we have batched experiments into 10 PsiCoverage files per tissue.


Individual steps with :bash:`majiq-build`
=========================================

The MAJIQ Builder pipeline chains together 3 different unique commands
from :bash:`majiq-build`:

.. command-output:: majiq-build
   :returncode: 1

You may want to use these commands directly rather than the pipeline because:

- Preprocessing GFF3 and SJ files separately for multiple analyses: this only
  needs to be done once per experiment/annotation and can be run as part of a
  standard RNA-seq pipeline after alignment.
- SJ files can be made in an embarassingly parallel manner on a
  cluster/distributed environment.
- Using the individual commands gives finer control as to where output files
  are saved.

It also gives access to :bash:`majiq-build combine`, which allows splicegraphs
from multiple analyses/different experiments to be combined/contrasted.
This is particularly useful for what we call a "two-pass build".


.. _quick-twopass:

Two-pass build
==============

The MAJIQ builder includes tools for contrasting and combining multiple
splicegraphs:

- :bash:`majiq-build combine` allows combining independent evidence from
  multiple splicegraphs into a single splicegraph.

    - This is roughly equivalent to running :bash:`majiq-build update` with the
      experiments from each build as independent build groups.
    - Largely negligible differences come from merging simplification calls on
      the component splicegraphs; if the differences are unacceptable, one can
      reset simplification calls with
      :bash:`majiq-build update --reset-simplify --simplify-only`.

- :bash:`majiq-build combine` allows treating novel junctions from some of
  these splicegraphs as known, highlighting junctions that were novel to
  specific experiments.
- :bash:`majiq psi-coverage`, with the :bash:`--ignore-from sg` flag, allows
  producing coverage for events that are unique to only one splicegraph
  (i.e. if it was structurally the same in the first build, ignore it).

    - This enables focusing on structurally novel events.
    - It can also prevent duplicate work with shared experiments/events which
      were quantified in previous builds that share the same events.

This functionality is of particular interest for our clinical analysis
pipelines for patients with suspected Mendelian disorders, where each
per-patient analysis shares a large group of controls.
The controls can be analyzed one time (first pass), with a second pass analysis
for each patient afterwards.
