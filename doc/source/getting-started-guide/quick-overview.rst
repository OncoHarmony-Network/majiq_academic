##############
Quick overview
##############

MAJIQ is a software package for defining and quantifying local splicing
variations (LSVs) from RNA-seq.
Conceptually, MAJIQ is divided into the following modules:

.. role:: bash(code)
   :language: bash

- MAJIQ Builder (:bash:`majiq-build`): Define splicegraphs and coverage over
  local splicing variations per input experiment.
- MOCCASIN (:bash:`majiq-moccasin`): Model observed and unobserved factors
  impact on input experiment coverage to obtain batch-corrected coverage.
- MAJIQ Quantifiers (:bash:`majiq-quantify`): Assess inclusion of spliced
  junctions and retained introns within/between groups of input experiments.
- MAJIQ Mendelian (working name) (:bash:`majiq-mendelian`): Summarize inclusion
  over control and case experiments to identify case-specific splicing outliers.

VOILA will be updated soon to enable visualization of spicing analysis with
this updated version of MAJIQ.
