

Output
======

In the /output/directory/ (specific by -d), the following files are
generated:

   `summary.tsv <#j38ua9j50u8t>`__

   `cassette.tsv <#drha9ovjrh8z>`__

   `tandem_cassette.tsv <#kix.i0mwd5sviwug>`__

   `alt5prime.tsv <#fh6j3l4uswjd>`__

   `alt3prime.tsv <#fh6j3l4uswjd>`__

   `alt3and5prime.tsv <#fh6j3l4uswjd>`__

   `alternate_first_exon.tsv <#y8aupght5kkv>`__

   `alternate_last_exon.tsv <#y8aupght5kkv>`__

   `p_alt5prime.tsv <#p4kscqdsxln2>`__

   `p_alt3prime.tsv <#p4kscqdsxln2>`__

   `mutually_exclusive.tsv <#ut40m3st4e53>`__

   `alternative_intron.tsv <#alternative-intron>`__

   `exitron.tsv <#q6hx5sn47zjf>`__

   `other.tsv <#other>`__

   constitutive.tsv (optional: see --keep-constitutive and
   --show-all-modules)

   `p_multi_gene_region.tsv <#i8llr21goq02>`__ (optional: see
   --multi-gene-regions)

For most users, we think the most important file is the summary file.
Each row in the file corresponds to a `Module <#splicegraph-modules>`__,
and the columns are as follows:

1.  Module: <Gene ID>_<N>

2.  Gene ID

3.  Gene Name

4.  Chr

5.  Strand

6.  LSV ID(s): <semi-colon separated LSV IDs>

7.  Cassette: <Count>

8.  Tandem Cassette: <Count>

9.  Alt 3: <Count>

10. Alt 5: <Count>

11. P_Alt 3: <Count>

12. P_Alt 5: <Count>

13. Alt 3 and Alt 5: <Count>

14. MXE: <Count>

15. Alternative Intron: <Count>

16. ALE: <Count>

17. AFE: <Count>

18. P_ALE: <Count>

19. P_AFE: <Count>

20. Orphan Junction: <Count>

21. Constitutive Junction: <Count> (Optional: see
       `--keep-constitutive <#bieyaxqd0clv>`__)

22. Constitutive Intron: <Count> (Optional: see
       `--keep-constitutive <#bieyaxqd0clv>`__)

23. Multi Exon Spanning: <Count>

24. Exitron: <Count>

25. Complex: <TRUE or FALSE>

..

   *Complex=False only if modules are comprised of a single event, OR
   meet specific*

   *criteria such as: (1) only comprised of 1 Intron Retention and 1
   P_Alt 5/3, or (2) only*

   *comprised of 1 Tandem Cassette, 1 Multi Exon Spanning, and any
   number of*

   *Constitutive Junctions*

26. Number of Events <Sum of Counts>

..

   *Sum of the event <Counts> across this module.*

27. Collapsed Event Name

*Summarizes the types and counts of events in this module.*

   *Tip: You can use Excel Pivot Table (Fields: ‘Module ID’ and
   ‘Collapsed Event Name’,*

   *Rows: Collapsed Event Name, Values: Count of Module ID) to tabulate
   of all the*

   *different types of modules in the summary.tsv file.*


