Installation of `mutpath`_
==============================

Requirements
---------------
You will need the following software to run `mutpath`_.

* `mutpath`_ is written in `Python`_. It has been tested on Mac OS X 10.6.8 with `Python`_ version 2.6.7, but it should probably work with other Python 2.6 and 2.7 versions as well. 

* `mutpath`_ uses `Python`_ to do some optional graphing. To do this graphing, you need to install `matplotlib`_. `mutpath`_ has been tested with `matplotlib`_ 1.0.0.

* Some of the scripts using `MUSCLE`_ to do sequence alignments. If you want to do these alignments, install `MUSCLE`_. `mutpath`_ has been tested with `MUSCLE`_ version 3.8.

* You will need to install `GraphViz`_ to visualize the mutational trajectories. `mutpath`_ has been tested with `GraphViz`_ version 2.30.1 on Mac OS X 10.6.8.

* You will need `BEAST`_ and `BEAGLE`_ in appropriate versions as detailed in the next section.


Installation
----------------

To install `mutpath`_, first download the source ZIP repository `on GitHub`_. After unzipping the file, run the following commands::

    cd mutpath
    python setup.py build
    python setup.py test
    python setup.py install

The last command might require you to use ``sudo`` if you do not have privileges to the default installation directory. Alternatively, install it locally with::

    python setup.py install --user

These commands install the `Python`_ modules and also install several scripts, which provide the most convenient high-level interface into the `mutpath`_ package. These scripts should be installed to be executable; they are present in the ``scripts`` subdirectory. The source code for the package is the in the ``src`` directory. There is also an ``examples`` subdirectory that contains some examples. Looking at these examples may be helpful for understanding the typical workflow. 

.. include:: weblinks.txt
