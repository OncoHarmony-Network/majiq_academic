.. currentmodule:: new_majiq

#############
API reference
#############

This page provides an auto-generated summary of MAJIQ's API. For more details
and examples, refer to the relevant chapters in the main part of the
documentation.


Primary classes/functions
=========================

.. autosummary::
   :toctree: generated/

   rng_resize
   rng_seed
   SJExperiment
   SpliceGraph
   Events
   SpliceGraphReads
   PsiCoverage
   DPsiPrior
   DeltaPsi
   Heterogen
   PsiControlsSummary
   PsiOutliers.from_psicov


SpliceGraph
===========

SpliceGraph I/O
---------------

.. autosummary::
   :toctree: generated/

   SpliceGraph.from_gff3
   SpliceGraph.from_zarr
   SpliceGraph.to_zarr
   SpliceGraph.to_sqlite


SpliceGraph components
----------------------

.. autosummary::
   :toctree: generated/

   SpliceGraph.contigs
   SpliceGraph.genes
   SpliceGraph.exons
   SpliceGraph.introns
   SpliceGraph.junctions
   SpliceGraph.exon_connections
   Contigs
   Genes
   Exons
   GeneIntrons
   GeneJunctions
   ExonConnections
