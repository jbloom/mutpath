Typical workflow
=================================================================
``mutpath`` installs a number of scripts that perform the main operations of the package. A mutational path can be built by using these scripts along with `BEAST`_. By default, the ``mutpath`` scripts are installed so that they can be run directory from the command line, taking a user-specified and created input text file.

The scripts themselves are documented individually in the following sections. This section gives an overall workflow of how these scripts can be used in conjunction with `BEAST`_.

Create the `BEAST`_ XML file
-------------------------------
First, use ``mutpath_parse_to_beastxml.py`` to create an input `BEAST`_ XML file from the FASTA sequence files. You can do this by manually creating the text input file and then running::

    mutpath_parse_to_beastxml.py infile.txt

where ``infile.txt`` is a text file that you have created that specifies the input variables. Note that this script will create files that provide information about sequences that are outliers in terms of their divergence. If you think that these might be mis-annotated sequences that should not be part of the tree, specify them for exclusion using the options provided in ``infile.txt`` and then run the script again. The final output of interest will be a `BEAST`_ XML file.

Run `BEAST`_
-----------------
Next, you will want to run `BEAST`_ on the XML file that you have created. Note that you will need to use `BEAGLE`_. Below is an example command if you have installed both `BEAST`_ and `BEAGLE`_ locally, and have an XML file named ``infile.xml``::

    java -Xmx2024m -Xms2024m -Djava.library.path=/home/jbloom/BEAGLE_libs/lib -cp ~/BEAST/build/dist/beast.jar dr.app.beast.BeastMain -beagle infile.xml

This will create a .trees and .log file with the prefix specified in ``infile.xml``. You may want to run several replicates of `BEAST`_ in different directories to create more MCMC data to analyze. Use `Tracer`_ to examine the ``.log`` files to see if you run enough steps of MCMC to give good effective sample sizes.


Compact the ``.trees`` file
--------------------------------
Running `BEAST`_ using the Markov Jumps mutation mapping will create a very large output ``.trees`` file. It may be helpful (although it is optional) to compact this file by removing unnecessary information. You can do this using ``mutpath_compact_trees.py`` by::

    mutpath_compact_trees.py treefile.trees

This will create a more compact file called ``treefile_compact.trees`` which you can use for all subsequent commands. You can delete the larger file using::

    rm treefile.trees

This command is not necessary, but will greatly reduce (typically by around 20-fold) the size of the ``.trees`` files.


Get the mutational paths
----------------------------------
After running `BEAST`_ and possibly compacing the mutational paths with ``mutpath_compact_trees.py``, you will have a large collection of trees in one or more ``.trees`` files. Use the script ``mutpath_get_paths.py`` to get the mutational paths from these trees. This script also merges all of the trees (without the mutation annotations) into a single ``.trees`` file that you can analyze with `TreeAnnotator`_ to make a maximum-clade credibility tree. The output of this script is a ``.txt`` file that contains the mutational paths between the specified starting and ending sequence. To do this, you manually create the text input file and then run::

    mutpath_get_paths.py infile.txt

The ``.txt`` file containing the mutational paths can then be analyzed by later scripts. 


Making and annotating the maximum clade credibility tree
-------------------------------------------------------------
After running ``mutpath_get_paths.py``, all of the various ``.trees`` files will have had their non-burnin trees combined into a single merged ``.trees`` file. Let's say we called that file ``merged_trees.trees``. You can then use `TreeAnnotator`_ to make a maximum clade credibility tree, which in turn can be visualized by `FigTree`_. To do this, make sure that you have `TreeAnnotator` installed such that the ``treeannotator`` executable is accessible, and then run::

    treeannotator merged_trees.trees maxcladecredibility.trees

to create the maximum clade credibility tree in ``maxcladecredibility.trees``. This tree can then be visualized with `FigTree`_. In order to format it for better viewing with `FigTree`_, you can use the script ``mutpath_annotate_tree.py``. Create an input file specifying how you want to label tip nodes in this tree, and then run::

    mutpath_annotate_tree.py infile.txt
    
This will create an annotated ``.tree`` file with the name that you specify in the input file. You can visualize this tree (and if you want further annotate it) using `FigTree`_.


Constructing the mutational trajectory
----------------------------------------
After running ``mutpath_get_paths.py`` to create the ``.txt`` file containing the mutational paths, you can then use ``mutpath_make_digraph.py`` to visualize and analyze the mutational trajectory. The mutational trajectory is the path from the starting sequence to the ending sequence through sequence space, and the probability of different paths in the trajectory is taken as proportional to the number of paths in the ``.txt`` mutational path file that take that portion of the trajectory. To create these paths, create an input file. Then run the script::

    mutpath_make_digraph.py infile.txt

This will generate a `DOT`_ file containing the mutational trajectory as a directed graph through protein sequence space that can be visualized in `GraphViz`_. You can also use `GraphViz`_ to save this trajectory in a PDF or other file format. In addition, this script will output information about the dates of mutations and the times for which nodes persisted.


.. _`BEAGLE`: http://beast.bio.ed.ac.uk/BEAGLE
.. _`BEAST`: http://beast.bio.ed.ac.uk/Main_Page
.. _`Tracer`: http://beast.bio.ed.ac.uk/Main_Page
.. _`TreeAnnotator` : http://beast.bio.ed.ac.uk/TreeAnnotator
.. _`on GitHub`: https://github.com/jbloom/mutpath
.. _`Jesse Bloom`: http://labs.fhcrc.org/bloom/
.. _`matplotlib`: http://matplotlib.org/
.. _`MUSCLE`: http://www.drive5.com/muscle/
.. _`FigTree` : http://tree.bio.ed.ac.uk/software/figtree/
.. _`GraphViz`: http://www.graphviz.org/
.. _`DOT` : http://www.graphviz.org/doc/info/lang.html
