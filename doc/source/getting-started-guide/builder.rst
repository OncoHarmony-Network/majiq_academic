.. _builder:


.. _MAJIQ Builder full:

MAJIQ Builder In Depth
----------------------

**General note**: in addition to the example below, you can find help composing your own commands + config by using the :ref:`command-builder`.

*Creating the config file*

Majiq needs to know some basic information about your study, annotation, and experiments before it can run. Provide this information
by creating a config file for majiq in simple text .INI format , using the template specified below:

::

    [info]
    bamdirs=/data/MGP/ERP000591/bam[,/data/MGP2/ERP000591/bam]
    sjdirs=/sj
    genome=mm10
    strandness=None[|forward|reverse]
    [experiments]
    Hippocampus=Hippocampus1,Hippocampus2
    Liver=Liver1,Liver2
    [optional]
    Hippocampus1=strandness:None,
    Liver2=strandness:reverse,

*Info:*

This is the study global information needed for the analsysis. The mandatory fields are:

- bamdirs: Comma separated list of paths where the bam files are located. If two files with the same names are in different paths, the first usable path is chosen.
- sjdirs: Comma separated list of paths where the sj files from previous builder runs are located. Used when running majiq with ``--incremental``
- genome: Genome assembly
- strandness=forward[\|reverse\|none]: Some of the RNASeq preparations are strand specific. This preparations can be reverse strand specific[reverse], forward strand specific [forward], or non-strand specific [none]. This parameter is optional, in which case None is used as default.

*Experiments:*

This section defines the experiments and replicates that are to be analyzed. Each line defines a condition and its name can be customized in the following way:

<group name>=<experiment file1>[,<experiment file2>]

where the experiment file is the sorted bam filename inside one of the bamdir directories (excluding extension .bam).
MAJIQ expects to find within the same directory bam index files for each experiment with the format <experiment file>.bam.bai.

- Multiple replicates within an experiment should be comma separated with no spaces.

*Optional:*

This section allows the user to specify changes in the global parameters specific to single experiments. The syntax goes as follows, <experiment file1>=<option1>:<value1>,...,<optionN>:<valueN>

The user only need to add the experiments that have different parameter value than the global parameters, using this section only when is needed.

(Currently only strandness has been adapted to be overwrite using this optional section)




