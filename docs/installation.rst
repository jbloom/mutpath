Installation of ``mutpath``
==============================

``mutpath`` is written in Python. It requires no external Python packages outside of the standard Python library, although it can optionally do some graphing that requires installation of `matplotlib`_. ``mutpath`` has been tested on Mac OS X 10.6.8 with Python 2.6.7 and `matplotlib`_ 1.0.0, although it should work with other similar versions as well and on other platforms as well.

Some of the ``mutpath`` scripts use the `MUSCLE`_ alignment program (version 3.8 has been tested with ``mutpath``), which must be downloaded and installed separately.

To install ``mutpath``, first download the source ZIP repository `on GitHub`_. After unzipping the file, run the following commands::

    cd mutpath
    python setup.py build
    python setup.py test
    python setup.py install

The last command might require you to use ``sudo`` if you do not have privileges to the default installation directory. Alternatively, install it locally with::

    python setup.py install --user

These commands install the Python modules and also install several scripts, which provide the most convenient high-level interface into the ``mutpath`` package. These scripts should be installed to be executable; they are present in the ``scripts`` subdirectory. The source code for the package is the in the ``src`` directory. There is also an ``examples`` subdirectory that contains some examples. Looking at these examples may be helpful for understanding the typical workflow. 


.. _`BEAGLE`: http://beast.bio.ed.ac.uk/BEAGLE
.. _`BEAST`: http://beast.bio.ed.ac.uk/Main_Page
.. _`on GitHub`: https://github.com/jbloom/mutpath
.. _`Jesse Bloom`: http://labs.fhcrc.org/bloom/
.. _`matplotlib`: http://matplotlib.org/
.. _`MUSCLE`: http://www.drive5.com/muscle/
