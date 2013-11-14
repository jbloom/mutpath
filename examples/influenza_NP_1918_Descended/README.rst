===========================================
influenza NP 1918 descended
===========================================

This analysis is found in the ``./examples/influenza_NP_1918_Descendend/`` subdirectory of the main `mutpath`_ repository `on GitHub`_.  In repository `on GitHub`_, the large `BEAST`_ output ``*.trees`` files are not included due to their very large file size.

This analysis involves influenza NPs descended directly from the 1918 virus.
These are
human seasonal H1N1 (1918-1957, 1977-2008),
human H2N2 (1957-1968), humanH3N2 (1968-2013), and North American swine lineages.

This analysis utilizes the `mapmuts`_ and `mutpath`_ packages. 

The analysis proceeded as follows:

    1) Construction of sequence sets from these lineages where outlier / mis-annotated sequences have been removed. These aligned sequences are found in *Combined_NP_proteins.fasta*.

    2) Construction of phylogenetic trees with mutations mapped to branches for the sequences in *Combined_NP_proteins.fasta* using `BEAST`_ starting with the input file *Combined_NP_proteins.xml*. The construction of trees and mutational paths using the `mutpath`_ package. Note that the actual `BEAST`_ output files are not included in the repository due to the very large file size.


Construction of sequence sets
------------------------------------
The first step was the construction of aligned coding sequence sets for the human and swine lineages. The ultimate set of aligned sequences are in *Combined_NP_nts.fasta* and *Combined_NP_proteins.fasta*. 

Separate analyses were done for the human and swine lineages. The sets are constructed by the scripts *get_human_seqs.py* and *get_swine_seqs.py* which construct the files *Human_NP_nts.fasta*, *Human_NP_proteins.fasta*, *Swine_NP_nts.fasta*, and *Swine_NP_proteins.fasta*. These two files were then concatenated to create *Combined_NP_nts.fasta* and *Combined_NP_proteins.fasta*. These scripts use the following input files:

    * *NPseqs.fasta* is the set of all unique full-length influenza A coding DNA sequences as downloaded from the `Influenza Virus Resource`_ on June-25-2013.

    * *Aichi68-NP.fasta* is the coding DNA sequence for A/Aichi/2/1968 (H3N2) NP as taken from reverse-genetics plasmid pHWAichi68-NP.

    * *JVI_82_8947_Anomalies.txt* is a list of the strain names for the sequences identified as anomalous (either frozen in time or recombinant) in Appendices 1 and 2 of `Krasnitz et al, 2008`.

    * *JDB_Anomalies.txt* is a list of strain names that I (Jesse D. Bloom, JDB) have found appear to be anomalous based on their tree positioning using analyses with `RAxML`_ and `Path-O-Gen`_ as described below.

Human lineages
~~~~~~~~~~~~~~~~~~~
The retained sequences are all human NPs that are descended directly from the 1918
viruses in the human lineage. This is H1N1 from 1918 to 1957, then
H2N2 from 1957 to 1968, and then H3N2 from 1968 to 2013. It also includes
human H1N1 sequences isolated between 1977 and 2008 (human seasonal H1N1), with
24 years are subtracted from the year based on the idea that these viruses 
appear to have been frozen since 1953 before reappearing as described 
in `dos Reis et al, 2009`_. 

Only sequences encoding unique proteins are retained, and a maximum of 5 sequences per year of any given subtype are retained to make the trees not too large. Anomalous sequences are removed, based on those specified by `Krasnitz et al, 2008`_, and also those that I find to be in strong violation of the molecular clock based on an analysis of non-date-stamped trees built with `RAxML`_ with `Path-O-Gen`_. Specifically, the analysis was done as follows:

    1) Extract the human NP sequences descended from the 1918 virus into *Human_NP_nts.fasta* and *Human_NP_proteins.fasta*::

        python get_human_seqs.py

    2) Build a `RAxML`_ tree (no date stamping) to allow visual inspection for possible outliers, and write this tree to the ``RAxML_output`` subdirectory::

        /Users/jbloom/standard-RAxML-master/raxmlHPC -w /Users/jbloom/mutpath/examples/influenza_NP_1918_Descended/RAxML_output -n Human_NP_nts -p 1 -m GTRCAT -s Human_NP_nts.fasta
        /Users/jbloom/standard-RAxML-master/raxmlHPC -w /Users/jbloom/mutpath/examples/influenza_NP_1918_Descended/RAxML_output -n Human_NP_proteins -p 1 -m PROTCATJTT -s Human_NP_proteins.fasta

       The trees were then visually inspected using `Path-O-Gen`_, and clear outliers from the molecular clock were removed by adding them to *JDB_Anomalies.txt* and re-running the analysis.

Swine lineages
~~~~~~~~~~~~~~
The retained sequences are swine NPs that are descended directly from the 1918.
These are taken only from North American viruses (USA, Canada; viruses from Mexico not included since many seem anomalous). According
to `Brockwell-Staats et al, 2009`_, the North American viruses are predominantly the classical H1N1 lineage from 1918 to 1998, and from then on the NP is maintained even as some of them reassort other viral genes.

Only sequences that encode unique proteins are retained. A maximum of 5 sequences per year of any given subtype are retained to make the trees not too large. Anomalous sequences are removed, based on those specified by `Krasnitz et al, 2008`_, and also those that I find to be in strong violation of the molecular clock based on an analysis of non-date-stamped trees built with `RAxML`_ with `Path-O-Gen`_. Also, there are some viruses that are not descended from the 1918 among the swine North American viruses -- they are probably mixed with other lineages or come from avian -- these are also removed. Specifically, the analysis was done as follows:

    1) Extract the swine NP sequences descended from the 1918 virus into *Swine_NP_nts.fasta* and *Swine_NP_proteins.fasta*::

        python get_swine_seqs.py

    2) Build a `RAxML`_ tree (no date stamping) to allow visual inspection for possible outliers, and write this tree to the ``RAxML_output`` subdirectory::

        /Users/jbloom/standard-RAxML-master/raxmlHPC -w /Users/jbloom/mutpath/examples/influenza_NP_1918_Descended/RAxML_output -n Swine_NP_nts -p 1 -m GTRCAT -s Swine_NP_nts.fasta
        /Users/jbloom/standard-RAxML-master/raxmlHPC -w /Users/jbloom/mutpath/examples/influenza_NP_1918_Descended/RAxML_output -n Swine_NP_proteins -p 1 -m PROTCATJTT -s Swine_NP_proteins.fasta

       The trees were then visually inspected using `Path-O-Gen`_, and clear outliers from the molecular clock were removed by adding them to *JDB_Anomalies.txt* and re-running the analysis.

Combined human and swine lineages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The human and swine files were combined into *Combined_NP_nts.fasta* and *Combined_NP_proteins.fasta*::

    cat Human_NP_nts.fasta Swine_NP_nts.fasta > Combined_NP_nts.fasta
    cat Human_NP_proteins.fasta Swine_NP_proteins.fasta > Combined_NP_proteins.fasta

`RAxML`_ was then used to build a tree (no date stamping) with::

        /Users/jbloom/standard-RAxML-master/raxmlHPC -w /Users/jbloom/mutpath/examples/influenza_NP_1918_Descended/RAxML_output -n Combined_NP_nts -p 1 -m GTRCAT -s Combined_NP_nts.fasta
        /Users/jbloom/standard-RAxML-master/raxmlHPC -w /Users/jbloom/mutpath/examples/influenza_NP_1918_Descended/RAxML_output -n Combined_NP_proteins -p 1 -m PROTCATJTT -s Combined_NP_proteins.fasta

This tree was then visually analyzed using `Path-O-Gen`_ to confirm that it appears to be fairly clock-like given the date stamps for the tips.


Mapping the mutational trajectory
------------------------------------

Running `BEAST`_
~~~~~~~~~~~~~~~~~~~~~

The file *Combined_NP_proteins.xml* was constructed from the sequences in *Combined_NP_proteins.fasta* as a `BEAST`_ input file using a combination of `BEAUTI`_ and hand-annotation. This XML file specifies date-stamped sequences, a strict molecular clock, and a JTT model of substitution.

This file was then used as the input for four different runs of `BEAST`_ (version 1.8pre Prelease r5356) using the `BEAGLE`_ (revision 1093) library, which were performed in the subdirectories ``run1/``, ``run2/``, etc. These runs were performed on the FHCRC's rhino cluster using ``sbatch`` with the command::

    sbatch run.sbatch
    
where the contents of the ``run.sbatch`` file was as follows::

    #!/bin/sh
    #SBATCH
    #PBS -l walltime=480:00:00
    echo "Starting..."
    java -Xmx4048m -Xms4048m -Djava.library.path=/home/jbloom/BEAGLE_libs/lib -cp ~/BEAST/build/dist/beast.jar dr.app.beast.BeastMain -beagle Combined_NP_proteins.xml
    echo "Finished."

The identical command was executed in all four run directories.

Inspection of the ``.log`` files with `Tracer`_ indicated that the runs (each of 20 million steps with trees saved every 10,000 steps) appeared to have equilibrated after about 2.5 million steps (the first 250 saved trees). If these are removed as burn-in and the four runs are combined, the effective sample sizes seem adequate to suggest MCMC convergence.

Each of the ``.trees`` files were compacted::

    mutpath_compact_trees.py run1/Combined_NP_proteins.trees
    mutpath_compact_trees.py run2/Combined_NP_proteins.trees
    mutpath_compact_trees.py run3/Combined_NP_proteins.trees
    mutpath_compact_trees.py run4/Combined_NP_proteins.trees

This created the files ``run1/Combined_NP_proteins_compact.trees``, etc.

Note that these ``.trees`` files are not included in the `mutpath`_ repository on GitHub due to large file sizes.

Building the trajectories
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three trajectories were then built:

    1) For human H3N2, the trajectory from A/Aichi/2/1968 (H3N2) to A/Texas/JMM_49/2012 (H3N2).

    2) For human (seasonal) H1N1, the trajectory from A/Memphis/13/1978 (H1N1) to A/Taiwan/11526/2008 (H1N1).

    3) For swine, the trajectory from A/swine/Wisconsin/1/1957 (H1N1) to A/swine/Indiana/A00968365/2012 (H1N1).

The human H3N2 trajectory was built using the command::

    mutpath_get_paths.py get_paths_infile_human_H3N2.txt
    mutpath_make_digraph.py make_digraph_infile_human_H3N2.txt

where the contents of ``get_paths_infile_human_H3N2.txt`` are::

    # input file to mutpath_get_paths.py
    intreefiles run1/Combined_NP_proteins_compact.trees run2/Combined_NP_proteins_compact.trees run3/Combined_NP_proteins_compact.trees run4/Combined_NP_proteins_compact.trees 
    burnin 250
    mergedtreesfile merged_Combined_NP_proteins_compact.trees
    fastafile Combined_NP_proteins.fasta
    seqtype protein
    startseq 1968.50_COUNT1_STRAIN_A/Aichi/2/1968_HOST_Human_SUBTYPE_H3N2_DATE_1968.50
    endseq 2012.93_COUNT98_STRAIN_A/Texas/JMM_49/2012_HOST_Human_SUBTYPE_H3N2_DATE_2012.93
    mutpathsfile human_H3N2_mutpaths.txt

and the contents of *make_digraph_infile_human_H3N2.txt* are::

    # input file to mutpath_make_digraph.py
    mutpathfile human_H3N2_mutpaths.txt
    translateseqs False
    dotfile human_H3N2_trajectory.dot
    minweight 0.01
    labelcutoff 0.6
    nodenamefile None
    mutationdates human_H3N2_mutationdates
    lasttipdate 2012.93
    persistencefile human_H3N2_persistence.txt

The human H1N1 trajectory was built similarly, using::

    mutpath_get_paths.py get_paths_infile_human_H1N1.txt
    mutpath_make_digraph.py make_digraph_infile_human_H1N1.txt

where the input files *get_paths_infile_human_H1N1.txt* and *make_digraph_infile_human_H1N1.txt* are modified to specify the correct human H1N1 sequences and dates.

The swine trajectory was built using::

    mutpath_get_paths.py get_paths_infile_swine.txt
    mutpath_make_digraph.py make_digraph_infile_swine.txt
    
where the input files *get_paths_infile_swine.txt* and *make_digraph_infile_swine.txt* are modified to specify the correct swine sequences and dates.

The key output of these runs are the `DOT`_ files displaying the trajectories, which can be visualized using `GraphViz`_::

    human_H1N1_trajectory.dot 
    human_H3N2_trajectory.dot 
    swine_trajectory.dot

These `DOT`_ files were opened with `GraphViz`_ and used to save PDF and JPG files::

    human_H1N1_trajectory.pdf
    human_H3N2_trajectory.pdf 
    swine_trajectory.pdf
    human_H1N1_trajectory.jpg
    human_H3N2_trajectory.jpg 
    swine_trajectory.jpg

These images are shown below.

swine influenza mutational trajectory
***************************************

.. figure:: ../examples/influenza_NP_1918_Descended/swine_trajectory.jpg
   :align: center
   :alt: swine_trajectory.jpg
   :width: 45%

   Mutational trajectory for swine influenza.


Human H1N1 mutational trajectory
***********************************

.. figure:: ../examples/influenza_NP_1918_Descended/human_H1N1_trajectory.jpg
   :align: center
   :alt: human_H1N1_trajectory.jpg
   :width: 55%

   Mutational trajectory for human H1N1.

Human H3N2 mutational trajectory
***********************************

.. figure:: ../examples/influenza_NP_1918_Descended/human_H3N2_trajectory.jpg
   :align: center
   :alt: human_H3N2_trajectory.jpg
   :width: 55%

   Mutational trajectory for human H3N2. Note that the beginning of this trajectory is slightly different from that in `Gong et al, 2013`_ possibly because of the inclusion of additional sequences from H2N2 that contribute to the early part of the phylogenetic tree.





Building the maximum clade credibility tree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition, the ``mutpath_get_paths.py`` runs created the merged ``.trees`` file *merged_Combined_NP_proteins.fasta*, which was used to build the maximum clade credibility tree *maxcladecredibility.trees* using `TreeAnnotator`_ with the command::

    ~/BEASTv1.7.5/bin/treeannotator merged_Combined_NP_proteins_compact.trees maxcladecredibility.trees

This tree was then manually edited using ``mutpath_annotate_tree.py`` to create ``annotated_maxcladecredibility.trees`` by the command::

    mutpath_annotate_tree.py annotate_tree_infile.txt

The output file ``annotated_maxcladecredibility.trees`` was then opened in `FigTree`_ where it was saved to ``handannotated_maxcladecredibility.trees`` and annotated further by hand. The image was then saved using `FigTree`_ as ``handannotated_maxcladecredibility.pdf`` and then converted to a JPG with::

    convert -density 400 handannotated_maxcladecredibility.pdf handannotated_maxcladecredibility.jpg
    
This tree is shown below:

.. figure:: ../examples/influenza_NP_1918_Descended/handannotated_maxcladecredibility.jpg
   :align: center
   :alt: handannotated_maxcladecredibility.jpg
   :width: 85%

   Maximum clade credibility tree of NPs descended from 1918 virus. The swine trajectory is in blue, the human H1N1 in green, and the human H3N2 in red.


.. _`mapmuts`: https://github.com/jbloom/mapmuts
.. _`mutpath`: https://github.com/jbloom/mutpath
.. _`Influenza Virus Resource`: http://www.ncbi.nlm.nih.gov/genomes/FLU/FLU.html
.. _`EMBOSS needle`: http://emboss.sourceforge.net/download/
.. _`Krasnitz et al, 2008`: http://jvi.asm.org/content/82/17/8947.abstract
.. _`BEAST`: http://beast.bio.ed.ac.uk/Main_Page
.. _`dos Reis et al, 2009`: http://www.ncbi.nlm.nih.gov/pubmed/19787384
.. _`TreeAnnotator`: http://beast.bio.ed.ac.uk/TreeAnnotator
.. _`FigTree` : http://tree.bio.ed.ac.uk/software/figtree/
.. _`RAxML` : https://github.com/stamatak/standard-RAxML
.. _`Path-O-Gen` : http://tree.bio.ed.ac.uk/software/pathogen/
.. _`Brockwell-Staats et al, 2009` : http://www.ncbi.nlm.nih.gov/pubmed/19768134
.. _`BEAGLE`: http://beast.bio.ed.ac.uk/BEAGLE
.. _`Tracer`: http://beast.bio.ed.ac.uk/Main_Page
.. _`BEAUTI`: http://beast.bio.ed.ac.uk/BEAUti
.. _`GraphViz`: http://www.graphviz.org/
.. _`DOT` : http://www.graphviz.org/doc/info/lang.html
.. _`Gong et al, 2013`: http://elife.elifesciences.org/content/2/e00631
.. _`on GitHub`: https://github.com/jbloom/mutpath
