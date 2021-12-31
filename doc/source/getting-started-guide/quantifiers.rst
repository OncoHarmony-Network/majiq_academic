.. _quantifiers:

#######################
Splicing quantification
#######################

.. role:: bash(code)
   :language: bash

MAJIQ has extensive functionality to quantify and assess changes in splicing.
This functionality starts with the PsiCoverage files produced by the MAJIQ
builder.
During the :bash:`majiq-build psi-coverage` step, MAJIQ identifies, per
experiment, which LSVs meet criteria for "quantifiability" by considering
thresholds on the minimum number of total reads and nonzero read positions.
Then, for quantifiable LSVs, MAJIQ infers posterior distributions of the
**P**\ ercent **S**\ pliced **I**\ n (PSI) for raw and bootstrap replicates of
LSV coverage as described in [Vaquero2016]_, [VaqueroAicherJewellGazzara2021]_.
