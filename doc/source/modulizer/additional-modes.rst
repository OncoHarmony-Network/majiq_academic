
Additional Arguments/Modes
==========================

Decomplexify stage arguments
----------------------------

+------------------------------------------------------------+
| | --decomplexify-reads-threshold <positive integer read #> |
| | Filter out junctions where the number of reads is        |
| | below a certain value (integer). If multiple input       |
| | files are used, only the biggest number of reads is      |
| | used. 0 means off.                                       |
|                                                            |
| Default: 1                                                 |
+------------------------------------------------------------+

+-----------------------------------------------------------+
| | --decomplexify-psi-threshold <float 0-1>                |
| | Filter out junctions where PSI is below a certain       |
| | value (between 0.0 and 1.0). If multiple input files    |
| | are used, only the highest PSI value is used. 0 is off. |
|                                                           |
| Default: 0.05                                             |
+-----------------------------------------------------------+

+----------------------------------------------------------+
| | --decomplexify-deltapsi-threshold <float 0-1>          |
| | Filter out junctions where abs(E(dPSI)) is below a     |
| | certain value (between 0.0 and 1.0). If multiple input |
| | files are used, only the biggest difference (dPSI)     |
| | value is used. 0 is off.                               |
|                                                          |
| Default: 0                                               |
+----------------------------------------------------------+

Definitions for changing and non-changing
-----------------------------------------

+----------------------------------------------------------+
| --changing-threshold <float 0-1>                         |
|                                                          |
| | \*This argument applies only to deltapsi.voila files\* |
| | dPSI threshold over which junctions are considered     |
|                                                          |
| changing, but not necessarily with high confidence.      |
|                                                          |
| However, this argument is linked to                      |
|                                                          |
| --between-group-dpsi-threshold. Changing one will        |
|                                                          |
| change the other, unless you manually set both.          |
|                                                          |
| Default: 0.2                                             |
+----------------------------------------------------------+

+----------------------------------------------------------+
| --probability-changing-threshold <float 0-1>             |
|                                                          |
| | \*This argument applies only to deltapsi.voila files\* |
| | Threshold for Prob(E(dPSI)) over which junctions are   |
|                                                          |
| considered confidently changing.                         |
|                                                          |
| Default: 0.95                                            |
+----------------------------------------------------------+

--non-changing-threshold <float 0-1>

Default 0.05

Threshold for E(dPSI) under which junctions are considered non-changing

--probability-non-changing-threshold

Other Arguments
---------------

--show-all-modules

   By default, if a module lacks any LSVs, the module is discarded.
   Turning on this flag retains all modules, even if they lack LSVs.
   There are many genes and regions of genes where the read coverage was
   not sufficient for MAJIQ to quantify LSVs. These regions may contain
   bona-fide modules and alternative splicing events, but you may not
   want to consider these regions in your analysis, since our confidence
   in their quantifiability is low.

--keep-constitutive <positive integer read #>

   By our definition, modules occur when there are 2+ ways to create
   alternative isoforms. However, we think it will be useful for users
   to know where in a gene there are \*no\* alternative paths. Turning
   on this flag updates the definition of a Module to include splice
   junctions that constitutively splice two exons. Also, this flag
   categorizes an additional type of event in modules: `Constitutive
   Junctions <#n49xc0jzg0ym>`__. A given Constitutive Junction must have
   at least <# reads> in order to be categorized as constitutive. Unless
   defined by user, the default <# reads> is 1.

--putative-multi-gene-regions

Alternative mode for the voila Categorizer: identify where there are
breaks in gene splicegraphs

likely due to variable coverage across the gene, but possibly due to
multiple non-overlapping

transcripts derving from the same locus. See: `Putative Multi-Gene
Regions <#i8llr21goq02>`__.
