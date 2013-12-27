.. mutpath documentation master file, created by
   sphinx-quickstart on Thu Feb 21 15:37:33 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================================
`mutpath`_ overview
===================================

`mutpath`_ is a Python package for identifying mutational paths through sequence space generated using the MarkovJumps feature of `BEAST`_. 

This package was created by `Jesse Bloom`_. It is based on the approach used to generate the mutation trajectories in `Gong et al, 2013`_. The mutation-mapping `BEAST`_ itself was implemented by `Marc Suchard`_. It can be used to create mutational paths similar to that in `Figure 2A`_ of `Gong et al, 2013`_.

`mutpath`_ source code is available `on GitHub`_. You can download a ZIP file of the repository from this site by clicking on the ``Download ZIP`` button on the right side.

`mutpath`_ documentation is available `on GitHub Pages`_.

If you use `mutpath`_, please cite `Gong et al, 2013`_ as the reference describing the approach.

This package contains a number of scripts that interface with `BEAST`_. These scripts generate mutational trajectories through sequence space and other associated data from `BEAST`_ output. Briefly, `BEAST`_ is used to sample phylogenetic trees with mutations mapped to branches from the posterior distribution. These trees make it possible to calculate the posterior probability for different trajectories through sequence space between two known sequences.

Some example analyses are included in the ``./examples/`` subdirectory of the main package repository.



.. toctree::
   :maxdepth: 2

   installation
   beastbeagle
   workflow
   scripts
   examples
   pythonapi
   acknowledgements

.. include:: weblinks.txt
