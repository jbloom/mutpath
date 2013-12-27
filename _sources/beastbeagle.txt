BEAST and BEAGLE
=====================
`mutpath`_ is designed to interface with the `BEAST`_ phylogenetics package, so you will obviously want to install `BEAST`_. The actual mutational path mapping requires the Markov Jumps feature of `BEAST`_. This feature is not available in the current stable release of `BEAST`_ (this may change when version 1.8 is released), so you need to obtain the developer's version of `BEAST`_ (or perhaps version 1.8, I have not tested that). Use of the Markov Jumps feature also requires using the `BEAGLE`_ library, so you need to install that as well.

Here are the commands that I used to obtain `BEAST`_ (v1.8.0pre Prelease r5356) and `BEAGLE`_ (revision 1093) on my computer. For `BEAST`_::

    svn co http://beast-mcmc.googlecode.com/svn/trunk BEAST
    cd BEAST
    ant
    java -jar build/dist/beast.jar

And for `BEAGLE`_::

    svn co http://beagle-lib.googlecode.com/svn/trunk BEAGLE
    cd BEAGLE/
    ./autogen.sh
    ./configure
    make
    sudo make install
    make check

If you just want to install BEAGLE locally in some directory such as ``BEAGLE_libs``, use::

    ./configure --prefix ~/BEAGLE_libs

and then you can install without the ``sudo`` command.

.. include:: weblinks.txt
