=========================================
mutpath_annotate_tree.py 
=========================================

Script for annotating a `BEAST`_ ``.trees`` file (such as a maximum clade credibility tree from `TreeAnnotator`_) so that it is better formatted for visualization by `FigTree`_. 

To run this script, create an input file with the format described below. Then run this script by typing the script name followed by the name of the input file::

    mutpath_annotate_tree.py infile.txt

Or if the script is not executable on your platform::
    
    python mutpath_annotate_tree.py infile.txt

Motivation and details
-----------------------
`TreeAnnotator`_ can be used to create a ``.trees`` file that can be visualized with `FigTree`_. However, the visualization will typically not look very good. Although it can be manually edited within `FigTree`_, this script automates some of the annotations. Specfically, it allows you to show a set of user-specified names for just a subset of the branch tips. It also orders the tree nodes in ascending order.

Below is an example input file. After running the script on this input file, the new *outtreefile* file will be better formatted for opening by `FigTree`_. Specifically, tip node names will only be shown for tips listed in the input file.


Example input file
---------------------
The input file is a text file. Empty lines or lines beginning with # are ignored. Otherwise, the first line should have the key *intreefile* and then give the ``.trees`` file that we are annotated. The next line should have the key *outtreefile* and then give the name of the new output ``.trees`` file with the annotations. If *outtreefile* already exists, it is overwritten. Remaining lines should then give tip names as they are specified in *intreefile* followed by a space and the name that is given to that tip in the annotated tree. If a tip name is not listed, it is not annotated in the tree. Typically, if you are working with a large tree then you may want to annotate only a small subset of the tips. Here is an example input file::

    # input file for mutpath_annotate_tree.py
    intreefile prot_maxcladecredibility.trees
    outtreefile annotated_maxcladecredibility.trees
    A/Aichi/2/1968_1968.50 Aichi/1968
    A/Brisbane/10/2007_2007.10 Brisbane/2007

After running the script on this input file, two tips will be annotated as Aichi/1968 and Brisbane/2007, and the others will not be annotated.

.. include:: weblinks.txt
