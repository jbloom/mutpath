==========================================
mutpath_get_paths.py 
==========================================

Script for parsing mutational paths from ``.trees`` files created using `BEAST`_ with the Markov Jumps mutation-mapping feature. 

The ``.trees`` files contain trees with mutations mapped to branches as created by `BEAST`_. With specification of a starting and ending sequence (*startseq* and *endseq*), each tree represents a specific mutational path through sequence space from the starting to the ending sequence. For each tree (possibly excluding some initial trees as specified by *burnin* below), the mutational path is extracted and written to the text file *mutpathsfile*. In addition, all of the non *burnin* trees are also written without the branch mutation annotations to the new ``.trees`` file *mergedtreesfile* that can later be analyzed by the `TreeAnnotator`_.

The input ``.trees`` files can either be those directly produced by `BEAST`_, or can have been pre-processed by ``mutpath_compact_trees.py``. However, if they have not been pre-processed by ``mutpath_compact_trees.py`` then running this script may consume a very large amount of memory.

To run this script, first create a text infile of the format described
below in this documentation. If you name this input file ``infile.txt``, then you would 
run the script with the single
input file as an argument::

    mutpath_get_paths.py infile.txt

If the script is not executable on your platform, you can use::

    python mutpath_get_paths.py infile.txt


Input file format
-------------------
The input file should contain a series of keys (strings without spaces) that are followed by their values. Lines that begin with # are considered comments and are ignored, as are blank lines. The file should contain the following keys:

* *intreefiles* should specify a list of strings giving the names of one or more ``.trees`` files. These file names cannot contain spaces. These should be the ``.trees`` files created by running `BEAST`_ with the Markov Jumps mutation mapping. The ``.trees`` files can optionally have been compacted using ``mutpath_compact_trees.py``. Note that all of these ``.trees`` files must be for the same input files, and so should be identical up to the beginning of the tree entries.

* *burnin* should specify how many trees from each file in *intreefiles* are removed as burn-in (i.e. equilibration). Note that this is the number of trees that are removed, not the number of MCMC steps. So for example, if you ran an MCMC of 20,000,000 steps sampling at every 10,000 steps, there are 2,000 trees in each file. So t remove the first 10%, you would set *burnin* to 200.

* *mergedtreesfile* should specify the name of the new ``.trees`` file we create that contains all of the trees specified in *intreefiles* after removing the *burnin* trees from each file, and after removing the branch annotations. This file may be rather large if you have not first run ``mutpath_compact_trees.py`` on the *intreefiles*. If *mergedtreefiles* already exists, it is overwritten.

* *fastafile* should specify a FASTA file with headers equal to the sequence names used to create the `BEAST`_ input files used to generate the *intreefiles*, and then giving the full sequences. This file is used to annotate the sequences in the *intreesfiles*. In other words, when looking at the *intreefiles* taxa labels, each taxa label must correspond to a header in *fastafile*, and then have the sequence provided.

* *startseq* should specify the name of the starting sequence for the mutational path. This should be the same sequence name used in `BEAST`_ input file, and so should be present in the taxa labels in the ``.trees`` file. 

* *endseq* should specify the name of the ending sequence for the mutational path. Like *startseq*, it should be the same sequence name used in the taxa labels in the ``.trees`` file.

* *seqtype* should specify the sequence type: it can be "DNA" or "protein".

* *mutpathsfile* should specify the name of the new ``.txt`` file that is created and contains all of the mutational paths. The format is described below. If *mutpathsfile* already exists, it is overwritten.


Example input file
--------------------
Here is an example input file::

    # Input file to mutpath_get_paths.py
    #
    # List of six input files
    intreefiles prots_BEASTrun_1/prots_compact.trees prots_BEASTrun_2/prots_compact.trees prots_BEASTrun_3/prots_compact.trees prots_BEASTrun_4/prots_compact.trees prots_BEASTrun_5/prots_compact.trees prots_BEASTrun_6/prots_compact.trees 
    #
    # Remove the first 200 trees as burnin
    burnin 200
    #
    # File to write merged trees after removing burnin, branch annotations
    mergedtreesfile prot_trees_merged.trees
    #
    # FASTA file giving sequences for the strain names in intreefiles
    fastafile prots.fasta
    #
    # Sequence type: DNA or protein
    seqtype protein
    #
    # Name of starting and ending sequences for mutational path
    startseq A/Aichi/2/1968_1968.50
    endseq A/Brisbane/10/2007_2007.10
    #
    # File to write the mutational paths
    mutpathsfile prot_mutpaths.txt



Output files
----------------
Running ``mutpath_get_paths.py`` creates two output files, *mergedtreesfile* and *mutpathsfile*. 

* *mergedtreesfile* is a new ``.trees`` file that contains all of the trees in the *intreefiles*, after removing the first *burnin* trees from each file. It also contains all of the preface and suffix information necessary to allow it to be read by `TreeAnnotator`_. The mutation annotations on the branches have been removed, which shrinks the file to a manageable size (particularly if you have already run ``mutpath_compact_trees.py`` on the *intreefiles* to remove other information). You may want to use *mergedtreefiles* to generate a maximum-clade credibility tree with `TreeAnnotator`_.

* *mutpathsfile* is a new ``.txt`` file that contains the mutational paths. It lists each of the paths numbered starting at 1. Within each path, the mutations are indicated with numbering starting at 1 for the first position in the sequence. The times for the mutations, the starting and ending strains, and the most recent common ancestor of these two strains, are also indicated. These times are measured in units before the most recent tip node (so the root node would have the largest value of time). Here is an example of a file with two paths::

    MUTPATH 1
    startstrain_name A/Aichi/2/1968_1968.50
    startstrain_seq ATGGCAATGGGCTAA 
    startstrain_time 42.5
    endstrain_name A/Brisbane/10/2007_2007.10
    endstrain_seq ATGACGATTGGATAA
    endstrain_time 3.9
    commonancestor_seq ATGGCGATGGGCTAA
    commonancestor_time 43.12713
    startstrain_to_commonancestor_path
    A6G : 42.713     
    commonancestor_to_endstrain_path
    G9T : 31.732
    G4A : 25.1343
    C12A : 10.134

    MUTPATH 2
    startstrain_name A/Aichi/2/1968_1968.50
    startstrain_seq ATGGCAATGGGCTAA 
    startstrain_time 42.5
    endstrain_name A/Brisbane/10/2007_2007.10
    endstrain_seq ATGACGATTGGATAA
    endstrain_time 3.9
    commonancestor_seq ATGGCGATGGGCTAA
    commonancestor_time 44.12713
    startstrain_to_commonancestor_path
    A6G : 42.113     
    G9T : 43.124
    commonancestor_to_endstrain_path
    G4A : 21.1343
    C5A : 19.531
    A5C : 19.402
    C12A : 9.134



.. include:: weblinks.txt
