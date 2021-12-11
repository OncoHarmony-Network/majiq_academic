.. currentmodule:: new_majiq

#############
API reference
#############

This page provides an auto-generated summary of MAJIQ's API. For more details
and examples, refer to the relevant chapters in the main part of the
documentation.


Random number generation
========================

.. autosummary::
   :toctree: generated/

   rng_seed
   rng_resize


Build API
=========

Classes
-------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   SJExperiment
   SJJunctionsBins
   SJIntronsBins
   SpliceGraph
   Contigs
   Genes
   Exons
   GeneIntrons
   GeneJunctions
   ExonConnections
   ExperimentThresholds
   GroupJunctionsGenerator
   PassedJunctionsGenerator
   GroupIntronsGenerator
   SimplifierGroup


Create a splicegraph from GFF3
------------------------------

.. autosummary::
   :toctree: generated/

   SpliceGraph.from_gff3


Save/load splicegraphs to zarr
------------------------------

.. autosummary::
   :toctree: generated/

   SpliceGraph.to_zarr
   SpliceGraph.from_zarr


Process BAMs for junction/intron coverage
-----------------------------------------

.. autosummary::
   :toctree: generated/

   SJExperiment.from_bam
   SJExperiment.to_zarr
   SJExperiment.from_zarr
   SJExperiment.introns
   SJExperiment.junctions


Update SpliceGraph using intron, junction coverage
--------------------------------------------------


Update junctions
~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   GeneJunctions.builder
   GeneJunctions.build_group
   GroupJunctionsGenerator.add_experiment
   PassedJunctionsGenerator.add_group
   PassedJunctionsGenerator.get_passed


Update exons
~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   Exons.infer_with_junctions


Update introns
~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   Exons.empty_introns
   Exons.potential_introns
   GeneIntrons.update_flags_from
   GeneIntrons.build_group
   GroupIntronsGenerator.add_experiment
   GroupIntronsGenerator.update_introns
   GeneIntrons.filter_passed


Update SpliceGraph
~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   SpliceGraph.from_components
   SpliceGraph.with_updated_exon_connections
   ExonConnections.create_connecting


Update simplifier flags
-----------------------

.. autosummary::
   :toctree: generated/

   GeneIntrons._simplify_all
   GeneJunctions._simplify_all
   GeneIntrons._unsimplify_all
   GeneJunctions._unsimplify_all
   ExonConnections.simplifier
   SimplifierGroup.add_experiment
   SimplifierGroup.update_connections


Combine multiple splicegraphs together
--------------------------------------

.. autosummary::
   :toctree: generated/

   GeneJunctions.load_dataset
   GeneJunctions.combine_datasets
   GeneJunctions.from_dataset_and_genes


Events API
==========

Classes
-------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   Events
   UniqueEventsMasks

Create/save events objects
--------------------------

.. autosummary::
   :toctree: generated/

   ExonConnections.lsvs
   ExonConnections.constitutive
   PsiCoverage.get_events
   PsiControlsSummary.get_events
   Events.to_zarr
   Events.from_zarr

Work with events objects
------------------------

.. autosummary::
   :toctree: generated/

   Events.unique_events_mask
   Events.exons
   Events.introns
   Events.junctions
   Events.df
   Events.ec_dataframe

Information on unique events
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   Events.e_idx
   Events.ref_exon_idx
   Events.event_type
   Events.ec_idx_start
   Events.ec_idx_end
   Events.connections_slice_for_event

Information on connections per event
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   Events.ec_idx
   Events.is_intron
   Events.connection_idx
   Events.connection_gene_idx
   Events.connection_start
   Events.connection_end
   Events.connection_denovo
   Events.connection_ref_exon_idx
   Events.connection_other_exon_idx


PsiCoverage API
===============

Classes
-------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   PsiCoverage

Create/save PsiCoverage
-----------------------

.. autosummary::
   :toctree: generated/

   PsiCoverage.from_sj_lsvs
   PsiCoverage.to_zarr
   PsiCoverage.to_zarr_slice_init
   PsiCoverage.to_zarr_slice
   PsiCoverage.from_zarr
   PsiCoverage.updated
   PsiCoverage.sum
   PsiCoverage.mask_events

Events/prefixes with coverage
-----------------------------

.. autosummary::
   :toctree: generated/

   PsiCoverage.num_connections
   PsiCoverage.get_events
   PsiCoverage.num_prefixes
   PsiCoverage.prefixes
   PsiCoverage.event_passed
   PsiCoverage.num_passed
   PsiCoverage.passed_min_experiments

Raw coverage/posteriors
-----------------------

.. autosummary::
   :toctree: generated/

   PsiCoverage.raw_total
   PsiCoverage.raw_coverage
   PsiCoverage.raw_alpha
   PsiCoverage.raw_beta
   PsiCoverage.raw_psi_mean
   PsiCoverage.raw_psi_std
   PsiCoverage.raw_psi_mean_population_median
   PsiCoverage.raw_psi_mean_population_quantile

Bootstrap coverage/posteriors
-----------------------------

.. autosummary::
   :toctree: generated/

   PsiCoverage.num_bootstraps
   PsiCoverage.bootstrap_total
   PsiCoverage.bootstrap_coverage
   PsiCoverage.bootstrap_alpha
   PsiCoverage.bootstrap_beta
   PsiCoverage.bootstrap_psi_mean
   PsiCoverage.bootstrap_psi_mean_legacy
   PsiCoverage.bootstrap_psi_std
   PsiCoverage.bootstrap_psi_mean_population_median
   PsiCoverage.bootstrap_psi_mean_population_quantile

Beta approximation to bootstrap mixture coverage/posteriors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   PsiCoverage.approximate_alpha
   PsiCoverage.approximate_beta
   PsiCoverage.approximate_quantile
   PsiCoverage.approximate_discretized_pmf


Quantifier API
==============

DeltaPsi (replicate PsiCoverage)
--------------------------------

.. autosummary::
   :toctree: generated/

   DPsiPrior
   DPsiPrior.empirical_update
   DeltaPsi
   DeltaPsi.dataset
   DeltaPsi.discrete_posterior_mean
   DeltaPsi.discrete_posterior_std
   DeltaPsi.probability_changing
   DeltaPsi.probability_nonchanging
   DeltaPsi.discrete_bootstrap_posterior_mean
   DeltaPsi.discrete_bootstrap_posterior_std
   DeltaPsi.bootstrap_probability_changing
   DeltaPsi.bootstrap_probability_nonchanging


Heterogen (independent PsiCoverage)
-----------------------------------

.. autosummary::
   :toctree: generated/

   Heterogen
   Heterogen.dataset
   Heterogen.raw_stats
   Heterogen.bootstrap_stats
   Heterogen.approximate_stats


CLIN (in development)
---------------------

Controls
~~~~~~~~

.. autosummary::
   :toctree: generated/

   PsiControlsSummary
   PsiControlsSummary.from_psicov
   PsiControlsSummary.from_zarr
   PsiControlsSummary.to_zarr
   PsiControlsSummary.q
   PsiControlsSummary.num_passed
   PsiControlsSummary.prefixes
   PsiControlsSummary.passed_min_experiments
   PsiControlsSummary.psi_median
   PsiControlsSummary.psi_quantile
   PsiControlsSummary.psi_range

Outliers
~~~~~~~~

.. autosummary::
   :toctree: generated/

   PsiOutliers
