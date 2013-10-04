==========================================
``mutpath_make_digraph.py`` script
==========================================

Script for building mutational trajectory directed graphs.

Using the mutational paths file written by ``mutpath_get_paths.py``, this script builds a directed graph through sequence space showing the mutational
trajectory. This path can be visualized using `GraphViz`_, which is a freely available software package.

The mutational trajectory shows the evolutionary path through sequence space as
a directed graph consisting of a series of nodes (sequences) and edges (mutations
connecting specific nodes). The trajectory moves from one known sequence to
another known sequence, passing through some number of uncertain intermediates.
The trajectory is built from some number of possible mutational paths
from the starting to ending sequence. In the trajectory, the weight of each node
corresponds to the fraction of paths that contain that sequence, while the
weight of each edge corresponds to the fraction of paths that contain that edge.
Note that if a path contains a node or edge more than once (which can happen
if there are mutational cycles), the node or edge is only considered to have
occurred once in that path for the purposes of assigning the weights.

To run this script, create an input file of the format described below. Then
run the script with that input file as an argument::

    mutpath_make_digraph.py infile.txt

If the script is not executable on your platform, you can use::

    python mutpath_make_digraph.py infile.txt


Input file format
-------------------
The input file should contain a series of keys (strings without spaces) that are followed by their values. Lines that begin with # are considered comments and are ignored, as are blank lines. The file should contain the following keys:

* *mutpathfile* should specify the name of a file containing mutational paths. This file should be in the format created by ``mutpath_get_paths.py``. This file is used to generate the trajectory.

* *translateseqs* specifies whether we translate nucleotide sequences into protein sequences and build a trajectory in protein sequence space. If you set this to False, then the sequence paths in *mutpathfile* are NOT translated -- this would be appropriate if either these paths are already in protein sequence space, or they are in nucleotide sequence space but you want to build the trajectory in nucleotide sequece space. However, you may want to construct the paths from the nucleotide sequences, but then build the trajectory in protein sequence space. In that case, the sequences in *mutpathfile* should specify coding nucleotide sequences, and you should set *translateseqs* to True. The trajectory object will then consist of protein sequence nodes and edges constructed by translating these nucleotide sequences. Note that no check is performed to make sure that the sequences translate well (i.e. without stop codons, gaps, etc), so be careful you really have nucleotide sequences that are in-frame coding.

* *dotfile* should specify the name of the output ``.dot`` file created by the script. This file is written in the `DOT`_ language, and can be visualized using `GraphViz`_. If this file already exists, it is overwritten.

* *minweight* specifies the minimum weight (probability of being on the trajectory) that a node or edge must possess to be displayed. You might want to set this to a value of 0.01 or 0.01 to avoid making a hugh *dotfile* due to lots of very low weight nodes and edges. If you set it to zero, all nodes and edges with any nonzero weight will be included.

* *labelcutoff* specifies the minimum weight that an edge must possess in order for it to be labeled with text showing the mutation for that edge. In addition, all nodes with weights >= *labelcutoff* have an incoming and outgoing edge that is also labeled. If there is not an actual labeled graph edge with a weight >= *labelcutoff* to such nodes, then a line spanning several nodes is drawn in another color to connect. Generally, you want *labelcutoff* > 0.5, and it is possible that the program will run into problems if this is not the case.

* *nodenamefile* is a ``*.fasta`` file that allows you to specify names for specific nodes in the trajectory. These names are written on each specified node that also has a weight >= *labelcutoff* (nodes in *nodenamefile* with weights < *labelcutoff* are not labeled). If you do not want to use this option, then set this parameter to None. Otherwise, it should be the name of an existing text file. Each header should specify the name given to a node, and the sequences are the sequences of that node. Note that if *translateseqs* is True then *nodenamefile* should specify the names for the protein sequences, not the gene sequences, as the trajectory itself will be in protein sequence space in that case. For example, the file below specifies node names for four sequences, two as string names and two as numbers::

    >Aichi/1968
    MASTVW 
    >Brisbane/2007
    MATWVG
    >1 
    MATTVW
    >2 
    MATWVG

* *mutationdates* specifies the prefix for the files we use to display dates for all mutations that occur in at least *labelcutoff* of the paths. For each such mutation, the time that it first occurs along each path is recorded. The script then outputs the posterior median and the 90% Bayesian credible interval (median centered) of the time of occurrence for each of these mutations taken over the set of paths in which that mutation occurs. In addition, it indicates the posterior probability that the mutation occurs on the forward path from the common ancestor to the ending sequence (the posterior probability that the mutation is on the reverse path from the starting sequence to the common ancestor is one minus this number). In *mutpathfile* itself, times are measured in units before the most recent tip node on the trees. The times here are measured as *lasttipdate* minus those times, so if *lasttipdate* is set to the date of the last tip node, then the times listed here are actually dates. A text file with the prefix *mutationdates* and the suffix ``.txt`` lists this information. If `matplotlib`_ is available, a PDF file with the prefix *mutationdates* and the suffix ``.pdf`` is also created to plot this information. These files are overwritten if they already exist. If *translateseqs* is *True*, then these are the dates for the amino-acid mutations, not the nucleotide mutations.

* *lasttipdate* is a number that should give the date of the last tip node of the tree. The times in the mutation paths are then measured in units before this date. This date is used to get the absolute value of the dates in *mutationdates*.

* *persistencefile* is the name of a file to which we write the times that nodes persist before the trajectory moves to the next nodes. These persistence times are written only for nodes that occur and are present in *nodenamefile* and which also appear in a fraction of at least *labelcutoff* of the paths, and the nodes are labeled by the names specified in *nodenamefile*. The persistence of a node is the time until the move to the next node in the trajectory. If *translateseqs* is *False*, then this is the time until the next mutation. If *translateseqs* is *True*, then this is the time until the next nonsynonymous mutation. In this last case, the time until the next mutation of any type is also recorded. These times are only for paths that contain that node, and are the times until the first mutation for the first occurrence of the node in each path. The times are written as the posterior median over those paths containing the node, along with the 90% Bayesian credible interval (median centered).

Example input file
--------------------
Here is an example input file::

    # input file for mutpath_make_digraph.py
    mutpathfile prot_mutpaths.txt
    translateseqs False
    dotfile prot_trajectory.dot
    minweight 0.01
    labelcutoff 0.6
    nodenamefile nodenames.fasta
    mutationdates prot_mutationdates
    lasttipdate 2010.69
    persistencefile persistence.txt


Output files
----------------
This script creates the following output files:

* *dotfile* is a `DOT`_ formatted graph that can be visualized by `GraphViz`_. In the visualization, the areas of nodes are proportional to their weight, and the width of edges is proportional to their weight. The trajectory moves from the starting to ending sequence, with the position of nodes along the direction of movement (vertical) defined as the Hamming distance between that node and the starting sequence minus the Hamming distance between that node and the ending sequence plus the Hamming Distance between the starting and ending sequences.

* *mutationdates*.txt is a text file giving the first dates of occurrences of mutations that occur in at least *labelcutoff* of the paths.

* *mutationdates*.pdf is created if `matplotlib`_ is available. It is a PDF file plotting the first dates of occurrences of mutations that occur in at least *labelcutoff* of the paths. 

* *persistencefile* is a text file giving the persistence (time to next node) for nodes specified in *nodenamesfile* and that occur in at least *labelcutoff* of the paths.



.. _`BEAGLE`: http://beast.bio.ed.ac.uk/BEAGLE
.. _`BEAST`: http://beast.bio.ed.ac.uk/Main_Page
.. _`TreeAnnotator`: http://beast.bio.ed.ac.uk/TreeAnnotator
.. _`on GitHub`: https://github.com/jbloom/mutpath
.. _`Jesse Bloom`: http://labs.fhcrc.org/bloom/
.. _`matplotlib`: http://matplotlib.org/
.. _`MUSCLE`: http://www.drive5.com/muscle/
.. _`on GitHub Pages`: http://jbloom.github.com/mutpath/
.. _`GraphViz`: http://www.graphviz.org/
.. _`DOT` : http://www.graphviz.org/doc/info/lang.html
.. _`matplotlib`: http://matplotlib.org/

