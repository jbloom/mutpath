=============================================================
influenza human H3N2 NP
=============================================================

This analysis is found in the ``./examples/influenza_human_H3N2_NP/`` subdirectory in the ``mutpath`` package which can be obtained `on GitHub`_. This is an example analysis using ``mutpath`` on the nucleoprotein (NP) from human H3N2 influenza to build a mutational trajectory at the protein level.

This example essentially recapitulates the analysis in `Gong et al, 2013`_, although it is not the exact code used for that paper and so may not be a completely precise match.

In repository `on GitHub`_, the large `BEAST`_ output ``*.trees`` files are not included due to their very large file size.

Input sequence data
---------------------
The beginning data files in this directory are as follows:

``NPhumanH3N2.fasta`` is a FASTA file of all full-length coding nucleotide sequences from human H3N2 influenza, excluding lab strains, as downloaded from the Influenza Virus Resource on May-19-2011.

``Aichi1968_NP.fasta`` and ``Brisbane2007_NP.fasta`` are FASTA files containing the coding sequences from the NP genes from the A/Aichi/2/1968 and A/Brisbane/10/2007 H3N2 strains, as encoded in the Bloom lab reverse genetics plasmids pHWAichi68-NP (#593) and pHWBR10-NP (#589).


Generation of `BEAST`_ input XML file
----------------------------------------
The `BEAST`_ XML files were generated as follows using the ``mutpath_parse_to_beastxml.py`` script on my personal computer.

Input file ``mutpath_parse_to_beastxml-prot-infile.txt`` was created manually.

A files to specify excluded outlier sequences was created and named ``excluded_prots.txt``.

The script was then run::

    mutpath_parse_to_beastxml.py mutpath_parse_to_beastxml-prot-infile.txt

This created the ``.txt`` and ``.pdf`` files with prefixes ``prot-divergence_gaps`` / ``prot-divergence_identities``. These files were used to identify likely outlier sequences, which were added to the ``excluded_prots.txt`` file, and the script was run again.

The final output files were::

    prots.fasta
    prots.nex
    prots.xml

The XML files is the `BEAST`_ input file for the next step.


Running `BEAST`_
------------------
`BEAST`_ was used to run the ``prots.xml`` file. The version of `BEAST`_ used was v1.8.0pre Prerelease r5356, using the `BEAGLE`_ (revision 1093) library. This was done on the FHCRC's ``rhino`` computing cluster.

Six separate directories were created called ``prots_BEASTrun_1``, ``prots_BEASTrun_2``, ... The ``prots.xml`` file was copied into ``prots_BEASTrun_*``.
BEAST was then run in each of these directories. The runs were initiated by ``sbatch`` as described below. Note that it specifies the `BEAGLE`_ library location based on the fact that these have been installed locally in ``/home/jbloom/BEAGLE_libs``.

Specifically, to run `BEAST`_, the following file was created with the name ``runcommand.sbatch``::

    #!/bin/sh
    #SBATCH
    #PBS -l walltime=480:00:00
    echo "Starting..."
    java -Xmx4048m -Xms4048m -Djava.library.path=/home/jbloom/BEAGLE_libs/lib -cp ~/BEAST/build/dist/beast.jar dr.app.beast.BeastMain -beagle prots.xml > screenlog.txt
    echo "Finished."

The file was run using::

    sbatch runcommand.sbatch

This created the `BEAST`_ ``.trees`` and ``.log`` files.

The output directories were then copied back to my local computer for analysis. However, these directories are not included in the repository `on GitHub`_ due to the large file size.

The ``.log`` files were examined using `Tracer`_ to check for MCMC convergence. 

Using a burn-in of the first 2 million steps (first 10% of steps), the effective sample size (ESS) for both the posterior and the root height exceed 200, suggesting good MCMC convergence.


Compacting the ``.trees`` files
--------------------------------
The ``prot.trees`` files created by running `BEAST`_ were then compacted to new files called ``prots_compact.trees`` using the ``mutpath_compact_trees.py`` script.

The commands (run inside each individual run directory were as follows)::

    % du -h prots.trees
    939M    prots.trees
    % mutpath_compact_trees.py prots.trees
    % du -h prots_compact.trees
    40M prots_compact.trees
    % rm prots.trees


Extracting the mutational paths
----------------------------------
To extract the mutational paths and merged trees from the compacted ``.trees`` file, the script ``mutpath_get_paths.py`` was used. The first 10% of each `BEAST`_ run was specified as burnin. The file ``mutpath_get_prot-paths_infile.txt`` was created manually.
The script was then run using the following command::

    mutpath_get_paths mutpath_get_prot-paths_infile.txt

to extract the protein mutational paths. This created an output file listing the mutational paths (``prot_mutpaths.txt``) as well as the merged trees without branch annotations (``prot_trees_merged.trees``). The merged ``.trees`` file is designed for building maximum-clade credibility trees with `TreeAnnotator`_.


Making and annotating the maximum clade credibility tree
------------------------------------------------------------
The maximum clade credibility tree was constructed from the ``prot_trees_merged.trees`` files built by ``mutpath_get_trees.py`` using `TreeAnnotator`_. The commands was::

    ~/BEASTv1.6.1/bin/treeannotator prot_trees_merged.trees prot_maxcladecredibility.trees

to created the maximum clade credibility tree ``prot_maxcladecredibility.trees``.

This files was then further annotated with ``mutpath_annotate_tree.py``. An input file for this script was created and named ``mutpath_annotate_prot-tree.txt``, and then the script was run with::

    mutpath_annotate_tree.py mutpath_annotate_prot-tree.txt

This created the file ``annotated_prot_maxcladecredibility.trees``. `FigTree`_ was used to visualize this file, manually recolor the branches red, and save the image ``annotated_prot_maxcladecredibility.pdf``. This PDF image was used to create ``annotated_prot_maxcladecredibility.jpg`` using the shell ``convert`` utility.

.. figure:: ../examples/influenza_human_H3N2_NP/annotated_prot_maxcladecredibility.jpg
   :width: 65 %
   :align: center
   :alt: maximum clade credibility tree

   The maximum clade credibility tree (``annotated_prot_maxcladecredibility.jpg``).


Making the mutational trajectory and dating the mutations
------------------------------------------------------------
The mutational trajectory was created with ``mutpath_make_digraph.py`` from the ``prot_mutpaths.txt`` mutational paths file created by ``mutpath_get_paths.py``. The input file ``mutpath_make_digraph-prot_infile.txt`` was created manually, as was the file ``nodenames.fasta``, which defines the names of high-confidence nodes as used in "Stability-mediated epistasis constrains the evolution of an influenza protein." The ``mutpath_make_digraph.py`` script was then run on these input files with::

    mutpath_make_digraph.py mutpath_make_digraph-prot_infile.txt

This created a mutational trajectory in protein sequence space.
The trajectory was written in the `DOT`_ language in the file ``prot_trajectory.dot``.
It was then visualized using ``GraphViz`` (version 2.30), which was also used to save the image files ``prot_trajectory.pdf`` and ``prot_trajectory.jpg``. The script also created the following additional output files which contain information about the mutation dates and node persistence times::

    prot_mutationdates.pdf
    prot_mutationdates.txt
    prot_nodepersistence.txt

The ``prot_mutationdates.pdf`` file was used to create ``prot_mutationdates.jpg`` using the shell ``convert`` utility.

.. figure:: ../examples/influenza_human_H3N2_NP/prot_trajectory.jpg
   :width: 50 %
   :align: center
   :alt: the mutational trajectory

   The mutational trajectory (file ``prot_trajectory.jpg``).

.. figure:: ../examples/influenza_human_H3N2_NP/prot_mutationdates.jpg
   :width: 40 %
   :align: center
   :alt: the mutation dates

   The mutation dates (file ``prot_mutationdates.jpg``).



.. _`BEAGLE`: http://beast.bio.ed.ac.uk/BEAGLE
.. _`BEAST`: http://beast.bio.ed.ac.uk/Main_Page
.. _`Tracer`: http://beast.bio.ed.ac.uk/Main_Page
.. _`TreeAnnotator` : http://beast.bio.ed.ac.uk/TreeAnnotator
.. _`FigTree` : http://tree.bio.ed.ac.uk/software/figtree/
.. _`on GitHub`: https://github.com/jbloom/mutpath
.. _`Jesse Bloom`: http://labs.fhcrc.org/bloom/
.. _`GraphViz`: http://www.graphviz.org/
.. _`DOT` : http://www.graphviz.org/doc/info/lang.html
.. _`Gong et al, 2013`: http://elife.elifesciences.org/content/2/e00631
