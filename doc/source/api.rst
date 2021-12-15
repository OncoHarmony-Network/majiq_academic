.. currentmodule:: new_majiq

#############
API reference
#############

This page provides an auto-generated summary of MAJIQ's API. For more details
and examples, refer to the relevant chapters in the main part of the
documentation.


Random number generation
========================

MAJIQ uses a pool of random number generators to handle multithreaded random
number generation, which is separate from numpy or dask random number
generation.
If more than a single thread is needed for a task involving random numbers, do
not forget to use :py:meth:`~new_majiq.rng_resize` to size the pool of RNGs to
match.

.. autosummary::
   :toctree: generated/

   rng_seed
   rng_resize


Build API
=========

MAJIQ builds a :class:`SpliceGraph` object from GFF3 and coverage from BAMs.
The splicegraph is used later to define :class:`Events`: for quantification.

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
   GeneJunctionsAccumulator


Create a splicegraph from GFF3
------------------------------

The first step in MAJIQ is to build a splicegraph from transcriptome
annotations (GFF3).

.. autosummary::
   :toctree: generated/

   SpliceGraph.from_gff3


Save/load splicegraphs to zarr
------------------------------

These splicegraphs are saved and loaded with the following commands.

.. autosummary::
   :toctree: generated/

   SpliceGraph.to_zarr
   SpliceGraph.from_zarr


Process BAMs for junction/intron coverage
-----------------------------------------

Updating the annotated splicegraph requires information about coverage from
RNA-seq experiments, which is represented by
:py:class:`~new_majiq.SJExperiment` objects.

.. autosummary::
   :toctree: generated/

   SJExperiment.from_bam
   SJExperiment.to_zarr
   SJExperiment.from_zarr
   SJExperiment.introns
   SJExperiment.junctions


Update SpliceGraph structure, passed flags
------------------------------------------

SpliceGraphs are generally updated in the following manner:

- Update :py:class:`~new_majiq.GeneJunctions` using either coverage or existing
  :py:class:`~new_majiq.GeneJunctions` objects.
- Update :py:class:`~new_majiq.Exons` to match updated
  :py:class:`~new_majiq.GeneJunctions`.
- Get potential introns between exons. Update them using old introns and/or
  coverage, then filter for those that passed.
- Create updated :py:class:`~new_majiq.ExonConnections`, then
  :py:class:`~new_majiq.SpliceGraph`.


Update junctions using coverage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`GeneJunctions` can be updated using coverage from :class:`SJExperiment`
objects.
This is done by:

- Creating a :class:`PassedJunctionsGenerator` (:meth:`GeneJunctions.builder`).
- Experiments are passed in as "build groups" where evidence for a junction
  must be found in some minimum number of experiments of at least one build
  group.
- Build groups are represented by :class:`GroupJunctionsGenerator` (created by
  :meth:`GeneJunctions.build_group`).
- Add :class:`SJJunctionsBins` from experiments in a build group using
  :meth:`GroupJunctionsGenerator.add_experiment`, then add build groups using
  :meth:`PassedJunctionsGenerator.add_group`.
- The updated :class:`GeneJunctions` is then returned by
  :meth:`PassedJunctionsGenerator.get_passed`.

.. autosummary::
   :toctree: generated/

   GeneJunctions.builder
   GeneJunctions.build_group
   GroupJunctionsGenerator.add_experiment
   PassedJunctionsGenerator.add_group
   PassedJunctionsGenerator.get_passed


Update junctions from other junctions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Updated :class:`GeneJunctions` can also be created by loading junctions from
previous splicegraphs and combining them using :class:`GeneJunctions`.
Note that they must all share the same :class:`Genes` object, which can be done
by setting `genes` argument to :meth:`GeneJunctions.from_zarr`.

.. autosummary::
   :toctree: generated/

   GeneJunctionsAccumulator.add
   GeneJunctionsAccumulator.accumulated


Update exons
~~~~~~~~~~~~

Update :py:class:`~new_majiq.Exons` to match updated
:py:class:`~new_majiq.GeneJunctions`.

.. autosummary::
   :toctree: generated/

   Exons.infer_with_junctions


Update introns
~~~~~~~~~~~~~~

Generally, updated introns are obtained by:

- Determine potential introns between exons (:meth:`Exons.potential_introns`).
- Update flags of potential introns using old introns
  (:meth:`GeneIntrons.update_flags_from`) and/or coverage in build groups
  (:meth:`GeneIntrons.build_group`).
- Filter potential introns to only those that passed build filters
  (:meth:`GeneIntrons.filter_passed`)

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

The updated splicegraph is made by:

- Creating updated :py:class:`~new_majiq.ExonConnections`
  (:py:meth:`~new_majiq.ExonConnections.create_connecting`).
- Creating updated :py:class:`~new_majiq.SpliceGraph`
  (:py:meth:`~new_majiq.SpliceGraph.with_updated_exon_connections`).

.. autosummary::
   :toctree: generated/

   SpliceGraph.from_components
   SpliceGraph.with_updated_exon_connections
   ExonConnections.create_connecting


Update simplifier flags
-----------------------

Simplifier flags allow excluding introns and junctions that pass reliability
thresholds (raw readrates/nonzero bins) but have negligible coverage relative
to the events they are a part of (PSI).
These flags are updated in place by creating a :class:`SimplifierGroup`, which
accumulates :class:`SJExperiment` objects per group and updates intron/junction
flags for a group using :meth:`SimplifierGroup.update_connections`.


.. autosummary::
   :toctree: generated/

   GeneIntrons._simplify_all
   GeneJunctions._simplify_all
   GeneIntrons._unsimplify_all
   GeneJunctions._unsimplify_all
   ExonConnections.simplifier
   SimplifierGroup.add_experiment
   SimplifierGroup.update_connections


Events API
==========

An event is defined by a reference exon and connections (junctions and/or
intron) that all start or end at the reference exon.
The :class:`Events` class represents a collection of these events as arrays
over events (e_idx) and connections per event (ec_idx).
The mapping from events and event connections is specified by offsets yielding
the start/end indexes of ec_idx for each event.
The events use indexes to refer back to the splicegraph/exon connections that
were used to create them.

:class:`UniqueEventsMasks` and :meth:`Events.unique_events_mask` allow
identification of events that are unique or shared between two :class:`Events`
objects.
This has use for analyses involving multiple splicegraphs derived from a common
splicegraph (e.g. a common set of controls).

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

:class:`PsiCoverage` describes coverage over :class:`Events` in one or more
independent "prefixes".
:class:`PsiCoverage` can be created over events for a single experiment using
:meth:`PsiCoverage.from_sj_lsvs` (prefix is determined by prefix of original
BAM file, which is where "prefix" name originates).
New :class:`PsiCoverage` files can be subsequently created by loading them
together or aggregating coverage over multiple prefixes.
Finally, :class:`PsiCoverage` provides attributes and functions which enable
lazy computation of PSI posterior quantities using xarray/Dask.


Classes
-------

.. autosummary::
   :toctree: generated/
   :nosignatures:

   PsiCoverage

Create/save PsiCoverage
-----------------------

Create PsiCoverage from SJ coverage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   PsiCoverage.from_sj_lsvs

Save PsiCoverage to zarr
~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   PsiCoverage.to_zarr
   PsiCoverage.to_zarr_slice_init
   PsiCoverage.to_zarr_slice

Load and update PsiCoverage
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

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
