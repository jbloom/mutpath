=========================================
``mutpath_compact_trees.py`` script
=========================================

Script compacting `BEAST`_ ``.trees`` files produced using the Markov Jumps mutation mapping feature.

Run this script by typing the script name followed by the name of the ``.trees`` file that you want to compact, for example::

    mutpath_compact_trees.py largefile.trees

Or if the script is not executable on your platform::

    python mutpath_compact_trees.py largefile.trees

These commands will create a more compact version of the file called ``largefile_compact.trees``.

Motivation and details
-----------------------

Running `BEAST`_ while mapping mutations to the branches tends to create very large output ``.trees`` files due to the large amount of information written. Most of the information is not actually necessary for reconstructing the mutational path. This script makes compact versions by removing information that is not needed to build mutational paths using the ``mutpath`` package. You can then store these compact ``.trees`` files more easily.

This script takes as input a ``.trees`` file created by running `BEAST`_ with Markov Jumps mutation mapping. You provide it with a single argument giving the name of the non-compacted ``.trees`` file. It creates a new compact version of the file in which the suffix ``.trees`` has been replaced by ``_compact.trees``. If a file with the name of this compacted file already exists, the script exits printing an error message. If the specified input file does not exist, the script exists with an error.

Once you have created the ``_compact.trees`` file, you should be able to delete the original ``.trees`` file and still do all of the mutational mapping using the more compact version. However, the script does not delete the non-compacted file itself -- you must do this manually.


The ``.trees`` file is compacted by doing the following to all of the lines containing trees:

1) All of the entries of the form *c\_allTransitions=4.5* are removed.

2) All of the entries of the form *&states="MATGA"* (giving sequence states) are removed.

3) All of the *history_all* notations are replaced with the more compact *h* in the lines giving the mutations, sites, and times.

4) All of the decimal numbers specifying the branch lengths are rounded to just four decimal places (so *3.872141428850743* becomes *3.8721*, for example). Although this slightly decreases the numerical accuracy, for reasonable branch lengths the difference should be unimportant.

Example
----------

Here is example usage of this script to shrink a file called ``prots.trees`` into the more compact ``prots_compact.trees``, and then delete the larger ``prots.trees`` file::

    % du -h prots.trees
    939M    prots.trees
    % mutpath_compact_trees.py prots.trees
    % du -h prots_compact.trees
    40M prots_compact.trees
    % rm prots.trees


.. _`BEAGLE`: http://beast.bio.ed.ac.uk/BEAGLE
.. _`BEAST`: http://beast.bio.ed.ac.uk/Main_Page
.. _`on GitHub`: https://github.com/jbloom/mutpath
.. _`Jesse Bloom`: http://labs.fhcrc.org/bloom/
.. _`matplotlib`: http://matplotlib.org/
.. _`MUSCLE`: http://www.drive5.com/muscle/

