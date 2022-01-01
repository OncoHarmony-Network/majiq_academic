.. _quantifiers:

#######################
Splicing quantification
#######################

.. role:: bash(code)
   :language: bash

MAJIQ has extensive functionality to quantify and assess changes in splicing.
This functionality starts with the :py:class:`~new_majiq.PsiCoverage` files
produced by :bash:`majiq psi-coverage`.
For the quantifiable LSVs that this step identifies, MAJIQ infers posterior
distributions of the
**P**\ ercent **S**\ pliced **I**\ n (PSI) for raw and bootstrap replicates of
LSV coverage as described in [Vaquero2016]_, [VaqueroAicherJewellGazzara2021]_.

MAJIQ provides three different quantification modes:

+ :bash:`majiq psi`: summarize inclusion from a single group of experiments as
  a posterior distribution over PSI.
  Input experiments are either quantified independently or as replicates
  (depending on :bash:`--min-experiments` flag).
+ :bash:`majiq deltapsi`: summarize differences in inclusion between two groups
  as a posterior distribution over a difference in PSI (dPSI).
  Input experiments are considered as replicates (one distribution on PSI per
  group, combined for one distribution on dPSI between groups).
+ :bash:`majiq heterogen`: test differences in PSI between samples from groups.
  Input experiments are considered independent samples (one distribution on PSI
  *per sample*).
  Test for differences between groups using 2-sample statistical tests on
  posterior means and samples generated from posterior distributions.

These quantifiers generally produce an output TSV (to stdout, alternatively to
path indicated by :bash:`--output-tsv`).
One should also generally add the :bash:`--splicegraph SG` flag to annotate
the quantifications with information about the LSVs being quantified.
PsiCoverage is specified for :bash:`majiq psi` as positional arguments.
For the other two quantifiers, specify group names with :bash:`-n NAME1 NAME2`
and PsiCoverage for each group with :bash:`-grp{1,2} PSICOV [PSICOV ...]`.

Additional details on required/optional parameters can be found by adding
:bash:`--help` to the subcommand of interest
(e.g. :bash:`majiq deltapsi --help`).

Alternatively, these quantifiers wrap the :ref:`PsiCoverage <api-psicoverage>`
and :ref:`quantifier <api-quantifier>` APIs (in development), which can be used
directly for more control over the analysis.
