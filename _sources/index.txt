.. mutpath documentation master file, created by
   sphinx-quickstart on Thu Feb 21 15:37:33 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================================
``mutpath`` overview
===================================

``mutpath`` is a Python package for identifying mutational paths through sequence space generated using the MarkovJumps feature of `BEAST`_. The package also contains scripts that might be useful for generating `BEAST`_ XML files for other purposes.

``mutpath`` source code is available `on GitHub`_. You can download a ZIP file of the repository from this site by clicking on the ``Download ZIP`` button on the right side.

``mutpath`` documentation is available `on GitHub Pages`_.

This package contains a number of scripts that interface with `BEAST`_. These scripts generate mutational trajectories through sequence space and other associated data from `BEAST`_ output. Briefly, `BEAST`_ is used to sample phylogenetic trees with mutations mapped to branches from the posterior distribution. These trees make it possible to calculate the posterior probability for different trajectories through sequence space between two known sequences.

Some example analyses are included in the ``./examples/`` subdirectory of the main package repository.

This package was created by `Jesse Bloom`_. It is based on the approach used to generate the mutation trajectories in `Gong et al, 2013`_. The mutation-mapping `BEAST`_ itself was implemented by `Marc Suchard`_.


.. toctree::
   :maxdepth: 1

   installation
   beastbeagle
   workflow
   mutpath_parse_to_beastxml
   mutpath_compact_trees
   mutpath_get_paths
   mutpath_annotate_tree
   mutpath_make_digraph
   example_influenza_human_H3N2_NP
   example_influenza_NP_1918_Descended
   pythonapi
   acknowledgements


.. _`BEAGLE`: http://beast.bio.ed.ac.uk/BEAGLE
.. _`BEAST`: http://beast.bio.ed.ac.uk/Main_Page
.. _`on GitHub`: https://github.com/jbloom/mutpath
.. _`Jesse Bloom`: http://labs.fhcrc.org/bloom/
.. _`matplotlib`: http://matplotlib.org/
.. _`MUSCLE`: http://www.drive5.com/muscle/
.. _`on GitHub Pages`: http://jbloom.github.com/mutpath/
.. _`Marc Suchard`: http://faculty.biomath.ucla.edu/msuchard/
.. _`Gong et al, 2013`: http://elife.elifesciences.org/content/2/e00631
