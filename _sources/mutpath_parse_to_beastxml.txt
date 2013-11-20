==========================================
``mutpath_parse_to_beast_xml.py`` script
==========================================

Script for parsing sequences and building BEAST input XML file. You do not have to use this script if you have another preferred method of creating the XML file.

To run this script, first create a text infile of the format described
below in this documentation. If you name this input file ``infile.txt``, then you would 
run the script with the single
input file as an argument::

    mutpath_parse_to_beastxml.py infile.txt

If the script is not executable on your platform, you can use::

    python mutpath_parse_to_beastxml.py infile.txt

Various variables must be specified in the script input file as described below.

This script goes through a FASTA file and parses out the
desired sequences, and then uses them to create a BEAST input XML file
for running BEAST. You can use it to run BEAST with or without the MarkovJumps
feature (which allows mutational path mapping but increases the output file
sizes). There are options for running BEAST using protein sequences (with the JTT matrix)
or nucleotide sequences (with the HKY model).
Various other options include coalescent trees with
a constant population size using Jeffrey's prior, and a strict clock with
a gamma prior on the substitution rate.

This script also helps you to identify possible outlier sequences (see the
*divergenceranking* option), and you can then manually specify these
outliers for exclusion use the *excludeseqs* option. This requires
running the script once to identify the outliers, and then running
it again with *excludeseqs* set to those sequences.

The sequences are provided in the FASTA input files specified by *seqfile*, *startseq*, and *endseq* and should have headers of the following format::

    >cds:AEB71984 A/Boston/99/2009 2009/04/27 NP

or::

    >cds:AEB71984 A/Boston/99/2009 2009-04-27 NP

or::

    >cds:BAA00478 A/Texas/1/1977 1977-- NP

or::

    >cds:BAA00478 A/Texas/1/1977 1977// NP

The strain names are extracted and the dates appended after an underscore::
    
    A/Boston/99/2009_2009.32

or::

    A/Texas/1/1977_1977.0.50

In the creation of the output file specified by *outfileprefix*, any square brackets *[* or *]* are converted to regular parentheses *(* or *)* because square brackets cause problems for `TreeAnnotator`_ in later analyses. Any spaces in the sequence names are converted to underscores.

You are asked to specify a starting sequence *startseq* which serves as
the reference sequence for alignment and also provides an initial sequence
for the mutational path. You are also asked to provide an end sequence *endseq*
for the end of the mutational path -- if you do not actually plan to build
such a path, just make the *endseq* the same as the *startseq* and it will be
removed if *purgeredundant* is True.

The overall workflow of the script is as follows:

1) Combine all the sequences in *startseq*, *endseq*, and *seqfile*; convert to upper case.

2) Translate into protein sequences (if *translateseqs* is True).

3) Purge redundant sequences or substrings of others (if *purgeredundant* is True).

4) Align all sequences individually to *startseq*, stripping any gaps relative to
   this sequence.

5) Output .fasta, .nex, and .xml files using *outfileprefix* for naming. The
   strain names are written in the format of the strain name followed by an underscore and the date. If two sequences would
   have identical names in this format, then they are named as strain_2_date, etc.

The output of this script is the parsed sequences in the following formats:

* FASTA

* Nexus

* BEAST XML file

Input file description
-------------------------

The input file format is a .txt file that specifies the following entries:

*seqfile*  is the name of a FASTA file containing the large set of sequences to be analyzed.
These sequences should have headers in the following format::

    >cds:AEB71984 A/Boston/99/2009 2009/04/27 NP

or::

    >cds:BAA00478 A/Texas/1/1977 1977// NP

or::

    >cds:AEB71984 A/Boston/99/2009 2009-04-27 NP

In particular, the second entry is presumed to give the strain name,
and the third entry must give the date of isolation, including at least
a year and a month and day if available. The dates are important,
as date stamped sequences are used in BEAST. If a date only includes the
year, the numerical date assigned is halfway through that year (July 1).
Note that **stop codons must be removed** from the sequences in *seqfile* prior to running this script.

*startseq* should give the name of a FASTA file containing a single sequence
entry which is the starting sequence for the mutational trajectory. The
header should have the same format as for *seqfile*. This is also the reference
sequence to which all other sequences are aligned.

*endseq* should give the name of a FASTA file containing a single sequence
entry which is the ending sequence for the mutational trajectory. The header
should have the same format as for *seqfile*. If you do not have an end sequence
(are not making a mutational path) then just set this to the same as *startseq*,
and the *endseq* will be purged if *purgeredundant* is True.

*musclepath* should give the path to a directory that contains the MUSCLE
alignment program in the form of an executable of the name ``muscle``.

*translateseqs* is a Boolean switch specifying whether or not we translate
the sequences. If we do translate them, then all of the subsequent operations
(alignment, etc) are done on the translated sequences. An exception is raised
if any sequence cannot be translated. Can be set to either True or False.

*purgeredundant* is a Boolean switch specifying that we remove any sequences
that are identical or substrings of others. The order they are removed is
arbitrary, except that the sequences specified by *startseq* and *endseq*
are always retained. Can be set to either True or False.

*divergenceranking* specifies where we output information about the sequence
divergence of each sequence relative to *startseq*. After aligning each
sequence to *startseq* (first translating it if *translateseqs* is True), the
pairwise identity of each sequence relative to *startseq* is computed among
non-gap positions, and the fraction of gaps is computed. These identities
are written to the text files *divergenceranking*\_identities.txt and *divergenceranking*\_gaps.txt. Note that the listing
here is generated after the *purgeredundant* and *excludeseqs* have been done.
Plots of cumulative distributions for the identities and gaps are also
written to files *divergenceranking*\_identities.pdf and *divergenceranking*\_gaps.pdf if `matplotlib`_ is available. You can also set *divergenceranking* to None, in which case no divergence ranking is done.

*excludeseqs* specifies the name of a text file listing the names of
sequences that we are manually specifying for exclusion from the analysis.
For example, you might have identified these as outliers or otherwise suspicious
sequences. These outliers are removed after purging of redundant sequences, but
before generating the listing given by *divergenceranking*. This file
should list the full header for each sequence to be excluded on a separate line,
with the headers in the same format they appear in *seqfile*, *endseq*, or *startseq*. It is OK if there is more text after the header, meaning that you can copy directly from the *divergenceranking* text files. 
You can assign *excludeseqs* to either an empty but existing
file or to None if you don't want to exclude anything. For example, here is a valid file::

    cds:ADX87378 A/Beijing/080302/2009 2009/08/03 NP    0.93574
    cds:ADX21092 A/Cambodia/NHRCC00010/2009 2009/10/15 NP   0.93712
    cds:AAK18005 A/Hong Kong/497/97 1997// NP   0.93775
    cds:ACX46212 A/Darwin/4/2005 2005/06/01 NP  0.93775


*outfileprefix* is the prefix that is assigned to the .fasta, .nex, and .xml
files created by this script.

*chainlength* specifies the length of the MCMC chain in the BEAST xml file.

*screenloginterval* specifies the frequency (number of steps) for which output
is specified to be written to the screen in the BEAST xml file

*fileloginterval* specifies the frequency (number of steps) for which output
is specified to be written to the .log and .trees file in the BEAST xml file.

*usemarkovjumps* specifies whether we use the Markov Jumps feature to map mutations
onto branches. Doing this allows reconstruction of the mutational path, but
also increases the output file size. Can be True or False.

*sitemodel* specifies the site model of evolution that is used. Currently accepted
values are JTT (for protein sequences) or HKY (for nucleotide sequences). Make
sure the value that you set here is consistent with your choice for *translateseqs*.

Example input file
-----------------------

Here is an example input file. In the input files, any lines beginning with a # character or that are empty are ignored. Otherwise entries should be specified by the name of the variable being specified followed by the specified value::

    # input file for mutpath_parse_to_beastxml.py
    # Builds a BEAST XML file from translated nucleotide sequences.
    seqfile NPhumanH3N2.fasta
    startseq Aichi1968_NP.fasta
    endseq Brisbane2007_NP.fasta
    # path to MUSCLE on my computer
    musclepath /Users/jbloom/muscle3.8/
    translateseqs True
    purgeredundant True
    divergenceranking divergence
    excludeseqs excluded_sequences.txt
    outfileprefix NPhumanH3N2_aligned
    chainlength 20000000
    screenloginterval 1000
    fileloginterval 10000
    usemarkovjumps True
    sitemodel JTT

.. _`matplotlib`: http://matplotlib.org/
.. _`TreeAnnotator` : http://beast.bio.ed.ac.uk/TreeAnnotator
