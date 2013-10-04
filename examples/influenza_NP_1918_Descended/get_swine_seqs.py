"""Gets and aligns sequences for phylogenetic tree of 1918 descended swine NPs.

This is the counterpart to ``get_human_seqs.py`` for the 1918-descendend
swine NPs.

This script utilizes the `EMBOSS needle`_ program for the alignments,
and requires the ``mapmuts`` program.
"""


import sys
import os
import re
import random
import datetime
import mapmuts.align
import mapmuts.sequtils


def DateToOrdinal(datestring, refyear=1968):
    """Converts a date string to an ordinal date.

    *datestring* is a date given by a string such as '2007/2/13' (for
    Feb-13-2007), or '2007/2//' if no day is specified, or
    '2007//' if no day or month is specified. The '/' characters can
    also be '-'.

    *refdate* is an integer year from the approximate timeframe we are examining
    which is used to anchor the datestring date on the assumption
    that each year has 365.25 days.

    The returned value is a number (decimal) giving the date. If no
    day is specified, the 15th (halfway through the month) is chosen.
    If no month or day is specified, July 1 (halfway through the
    year) is chosen.

    >>> print "%.2f" % DateToOrdinal('2007/4/27', 1968)
    2007.32

    >>> print "%.2f" % DateToOrdinal('2007/4/', 1968)
    2007.29

    >>> print "%.2f" % DateToOrdinal('2007//', 1968)
    2007.50

    >>> print "%.2f" % DateToOrdinal('2007-4-27', 1968)
    2007.32

    """
    if not isinstance(refyear, int):
        raise ValueError('refyear is not an integer')
    refdate = datetime.date(refyear, 1, 1).toordinal()
    try:
        if '/' in datestring:
            (year, month, day) = datestring.split('/')
        else:
            (year, month, day) = datestring.split('-')
    except ValueError:
        raise ValueError("Invalid datestring of: %s" % str(datestring))
    if year and month and day:
        (year, month, day) = (int(year), int(month), int(day))
        date = datetime.date(year, month, day)
    elif year and month:
        (year, month) = (int(year), int(month))
        date = datetime.date(year, month, 15)
    elif year:
        year = int(year)
        date = datetime.date(year, 7, 1)
    else:
        raise ValueError("Invalid datestring of: %s" % str(datestring))
    return (date.toordinal() - refdate) / 365.25 + refyear


def ParseHeader(head):
    """Parses an influenza sequence header for year, subtype, and host.

    Only works when the segment is NP.

    *head* is a string specifying a sequence header, in a format such as
    these examples::

        cds:CAA24268 Human A/Puerto Rico/8/1934 1934// H1N1 Puerto Rico NP
        cds:BAN57584 Avian A/muscovy duck/Vietnam/LBM330/2013 2013/01/07 H5N1 Viet Nam NP
        cds:AGN51218 Swine A/swine/Illinois/A01158960/2013 2013/02/08 H1N1 USA NP

    This function parses the header and returns a 4-tuple 
    *(strain, host, date, subtype, country)*. The date is extracted using
    *DateToOrdinal*.

    For example, the three example header above would return (if dates
    are rounded to one decimal for display):

        * *('A Puerto Rico/8/1934', 'Human', 1934.5, 'H1N1', 'Puerto Rico')*
        
        * *('A/muscovy duck/Vietnam/LBM330/2013', 'Avian', 2013.0, 'H5N1', 'Vietnam')*
        
        * *('A/swine/Illinois/A01158960/2013', 'Swine', 2013.1, 'H1N1', 'USA')*

    If the header cannot be parsed, returns *None*. This happens when:

        * The header contains the string *mixed*.

        * No 4-digit year is specified (i.e. the date is less than 1800)

        * The subtype is not fully known

        * No host is specified

    """
    if 'mixed' in head:
        return
    headmatch = re.compile('^\S+ (?P<host>[\w ]*) (?P<strain>A(\/[\w\. \-\'\(\)\?\+]*){1,6}) +(?P<date>(\d*|unknown)\/\d*\/\d*) (?P<subtype>[\w\,]+) (?P<country>[\w ]*) NP$')
    m = headmatch.search(head)
    if not m:
        raise ValueError("Failed to match header:\n%s" % head)
    host = m.group('host')
    if not host:
        return
    date = m.group('date')
    if 'unknown' in date or date == '//':
        return
    date = DateToOrdinal(date)
    if date < 1800:
        return
    subtypematch = re.compile('^H\d+N\d+$')
    subtype = m.group('subtype')
    if not subtypematch.search(subtype):
        return
    country = m.group('country')
    return (m.group('strain'), host, date, subtype, country)


def GetUnique(seqs):
    """Gets unique sequences from a set.

    *seqs* is a list of sequences as *(header, sequence)* 2-tuples.

    This method returns a new list *uniqueseqs* which is like *seqs*
    except that non-unique sequences have been removed. Specifically,
    all sequences that are identical to another sequence are removed.
    The first sequence is retained.
    """
    seq_d = {}
    for (head, seq) in seqs:
        if seq not in seq_d:
            seq_d[seq] = head
    uniqueseqs = [(head, seq) for (seq, head) in seq_d.iteritems()]
    return uniqueseqs


def PairwiseStatistics(s1, s2):
    """Returns alignment statistics for pairwise aligned sequences.

    *s1* and *s2* are strings of the same length corresponding to two
    aligned sequences.

    The returned variable is the 2-tuple *(f_ident, f_gaps)*. *f_ident*
    is the fraction of identities among the aligned non-gap positions.
    *f_gaps* is the fraction of position in the alignment that are gaps
    in one of the sequences.
    """
    assert len(s1) == len(s2)
    n_ident = n_gaps = 0
    for (x1, x2) in zip(s1, s2):
        if x1 == '-' or x2 == '-':
            n_gaps += 1
        elif x1 == x2:
            n_ident += 1
    return (n_ident / float(len(s1)), n_gaps / float(len(s1)))


def NeedleCDSandProtAlignments(refseq, seqs, needlecmd, tempfile='_alignments.temp'):
    """Uses EMBOSS needle to align sequences in *seqs* to *refseq*.

    The sequences in *seqs* and *refseq* should be coding DNA sequences.
    They must be translatable, but ambiguous nucleotides and truncation
    of incomplete sequences are allowed.

    *refseq* is a string giving the reference sequence.

    *seqs* is a list of *(header, sequence)* 2-tuples.

    *needlecmd* is the path to the EMBOSS needle program.

    *tempfile* is the name of a file temporarily created to store
    the alignments, and then deleted.

    Returns the 2-tuple *(prot_alignment, cds_alignment)*.

    *prot_alignment* contains each of the proteins encoded in 
    *seqs* as a *(header, protsequence)* 2-tuple, with the
    protein sequence aligned to that in *refseq* and with all
    gaps relative to *resfeq* stripped away.

    *cds_alignment* is an alignment of the coding DNA sequences
    in *prot_alignment*, with the nucleotide alignments done
    according to the protein alignments.
    """
    prots = []
    heads = []
    refprot = mapmuts.sequtils.Translate([('head', refseq)])[0][1]
    for (head, seq) in seqs:
        try:
            (head, prot) = mapmuts.sequtils.Translate([(head, seq)], readthrough_n=True, truncate_incomplete=True)[0]
            heads.append(head)
            prots.append(prot)
        except:
            sys.stderr.write("PROBLEM translating sequence %s" % head)
            raise
    try:
        mapmuts.align.Needle(refprot, prots, needlecmd, 'protein', tempfile)
        alignments = mapmuts.sequtils.ReadFASTA(tempfile)
    finally:
        if os.path.isfile(tempfile):
            os.remove(tempfile)
    assert len(alignments) == 2 * len(prots) == 2 * len(heads) == 2 * len(seqs)
    prot_alignment = []
    cds_alignment = []
    for i in range(len(prots)):
        prot = prots[i]
        head = heads[i]
        seq = seqs[i][1]
        assert seqs[i][0] == head
        (refa, prota) = (alignments[2 * i][1], alignments[2 * i + 1][1])
        assert len(refa) == len(prota)
        iref = iprot = 0
        alignedprot = []
        alignedcds = []
        for (aa_ref, aa_prot) in zip(refa, prota):
            assert (aa_ref == '-' or aa_ref == refprot[iref])
            assert (aa_prot == '-' or aa_prot == prot[iprot])
            if aa_ref == '-' and aa_prot != '-':
                iprot += 1
            elif aa_prot == '-' and aa_ref != '-':
                alignedprot.append(aa_prot)
                alignedcds.append('---')
                iref += 1
            elif aa_ref != '-' and aa_prot != '-':
                alignedprot.append(aa_prot)
                alignedcds.append(seq[3 * iprot : 3 * iprot + 3])
                iref += 1
                iprot += 1
            else:
                raise ValueError("Both prots in alignment have gap")
        alignedprot = ''.join(alignedprot)
        alignedcds = ''.join(alignedcds)
        assert alignedprot == mapmuts.sequtils.Translate([(head, alignedcds)], readthrough_n=True, truncate_incomplete=True, translate_gaps=True)[0][1]
        assert len(alignedprot) == len(refprot)
        prot_alignment.append((head, alignedprot))
        cds_alignment.append((head, alignedcds))
    assert len(prot_alignment) == len(cds_alignment)
    return (prot_alignment, cds_alignment)



def main():
    """Main body of script."""
    random.seed(1)

    # input / output files and script parameters
    hosts = frozenset(['Swine']) # only keep sequences from hosts in this set
    countries = frozenset(['USA', 'Canada'])
    nperyearhostsubtype = 5 # keep this many per year / host / subtype
    nperyearhost = 50 # keep this many per host per year regardless of subtype
    purgeambiguous = True # remove sequences with ambiguous nucleotides
    needlecmd = '/Users/jbloom/EMBOSS-6.5.7/emboss/needle' # path to EMBOSS needle
    refseqfile = 'Aichi68-NP.fasta' # sequence to which we align
    infile = 'NPseqs.fasta' # input file with all sequences
    outfile = 'Swine_NP_nts.fasta' # output file
    outprotfile = 'Swine_NP_proteins.fasta' # output file
    anomalies = frozenset([line.strip() for line in open('JVI_82_8947_Anomalies.txt').readlines() + open('JDB_Anomalies.txt').readlines()]) # strain names of anomalous sequences
    onlyunique = True # only keep unique sequences
    onlyuniqueproteins = True # only keep sequences encoding unique proteins
    subtypes_to_retain = { # Only keep subtypes in this dictionary.
        # for each subtype, the values are a dictionary keyed by every
        # year for which we keep the sequences, and with the values
        # being a number subtracted off the year for that subtype / year.
        'H1N1':dict([(year, 0) for year in range(1918, 2014)]),
        'H3N2':dict([(year, 0) for year in range(1998, 2014)]),
        'H1N2':dict([(year, 0) for year in range(1998, 2014)]),
    }


    # get and align the sequences
    if not os.path.isfile(needlecmd):
        raise IOError("Could not find EMBOSS needle of %s" % needlecmd)
    if not os.path.isfile(refseqfile):
        raise IOError("Could not find refseqfile of %s" % refseqfile)
    refseq = mapmuts.sequtils.ReadFASTA(refseqfile)
    if len(refseq) != 1:
        raise IOError("Failed to read exactly one sequence from refseqfile %s" % refseqfile)
    refprot = mapmuts.sequtils.Translate(refseq)[0][1]
    (refhead, refseq) = refseq[0]
    print "\nUsing a reference sequence of %s.\nThis sequence is %d nucleotides in length, and encodes a protein of %d residues." % (refhead, len(refseq), len(refprot))
    if not os.path.isfile(infile):
        raise IOError("Could not find infile %s" % infile)
    seqs = mapmuts.sequtils.ReadFASTA(infile)
    random.shuffle(seqs) # so they are not in order in case the download ordered them in some way
    print "\nRead %d sequences from %s." % (len(seqs), infile)
    print "\nPurging pre-defined anomalous sequences..."
    cleanseqs = []
    for (head, seq) in seqs:
        for anomaly in anomalies:
            if anomaly in head:
                break
        else:
            cleanseqs.append((head, seq))
    seqs = cleanseqs
    print "Retained %d sequences after removing anomalous ones." % len(seqs)
    if purgeambiguous:
        print "\nPurging any sequences with ambiguous nucleotide identities."
        m = re.compile('^[ATCGatcg]+$')
        cleanseqs = []
        for (head, seq) in seqs:
            if m.search(seq):
                cleanseqs.append((head, seq))
        seqs = cleanseqs
        print "Retained %d sequences after purging those with ambiguous nucleotides." % len(seqs)
    if onlyunique:
        print "\nGetting unique sequences..."
        seqs = GetUnique(seqs)
        print "Retained %d unique sequences." % len(seqs)
    if onlyuniqueproteins:
        print "\nGetting sequences encoding unique proteins.."
        prots = mapmuts.sequtils.Translate(seqs, readthrough_n=True, truncate_incomplete=True)
        uniqueprots = dict(GetUnique(prots))
        seqs = [(head, seq) for (head, seq) in seqs if head in uniqueprots]
        print "Retained %d sequences encoding unique proteins." % len(seqs)
    print "\nMaking sure the host / year / subtype can be parsed..."
    seqs = [(ParseHeader(head), seq) for (head, seq) in seqs if ParseHeader(head)]
    print "Retained %d sequences for which this information could be parsed." % len(seqs)
    print "\nRetaining sequences for which the host is %s..." % (', '.join(hosts))
    cleanseqs = []
    hostcounts = dict([(host, 0) for host in hosts])
    for tup in seqs:
        if tup[0][1] in hosts:
            hostcounts[tup[0][1]] += 1
            cleanseqs.append(tup)
    seqs = cleanseqs
    seqs = [tup for tup in seqs if tup[0][1] in hosts]
    print "Retained %d sequences from these hosts.\nCounts per host are:" % len(seqs)
    for (host, n) in hostcounts.iteritems():
        print "  %s: %d" % (host, n)
    print "\nRetaining only sequences for which the country is %s..." % (', '.join(countries))
    seqs = [tup for tup in seqs if tup[0][-1] in countries]
    print "Retained %d sequences from these countries." % len(seqs)
    print "\nNow retaining just %d sequences per host / subtype / year, and %d sequences per host / year regardless of subtype..." % (nperyearhostsubtype, nperyearhost)
    print "\nAlso only retaining sequences of the following subtypes in their specified date ranges: %s" % (', '.join(subtypes_to_retain.keys()))
    cleanseqs = []
    retainedcounts = {}
    hostretainedcounts = {}
    hostcounts = dict([(host, 0) for host in hosts])
    count = 1
    for ((strain, host, date, subtype, country), seq) in seqs:
        year = int(date)
        key = (host, year, subtype)
        if subtype in subtypes_to_retain:
            if year in subtypes_to_retain[subtype]:
                date -= subtypes_to_retain[subtype][year]
                year -= subtypes_to_retain[subtype][year]
            else:
                continue
        else:
            continue
        if key in retainedcounts:
            if retainedcounts[key] < nperyearhostsubtype:
                retainedcounts[key] += 1
            else:
                continue
        else:
            retainedcounts[key] = 1
        key2 = (host, year)
        if key2 in hostretainedcounts:
            if hostretainedcounts[key2] < nperyearhost:
                hostretainedcounts[key2] += 1
            else:
                continue
        else:
            hostretainedcounts[key2] = 1
        hostcounts[host] += 1
        cleanseqs.append(('%.2f_COUNT%d_STRAIN_%s_HOST_%s_SUBTYPE_%s_DATE_%.2f' % (date, count, strain, host, subtype, date), seq))
        count += 1
    seqs = cleanseqs
    print "Retained %d sequences overall.\nCounts per host are:" % len(seqs)
    for (host, n) in hostcounts.iteritems():
        print "  %s: %d" % (host, n)
#    print "Now adding in the reference sequence %s" % (refhead)
#    parsedhead = ParseHeader(refhead)
#    parsedhead = '%.2f_COUNT%d_STRAIN_%s_HOST_%s_SUBTYPE_%s_DATE_%.2f' % (parsedhead[2], 0, parsedhead[0], parsedhead[1], parsedhead[3], parsedhead[2])    
#    seqs = [(parsedhead, refseq)] + seqs
    print "\nNow translating and aligning..." 
    (prot_alignments, cds_alignments) = NeedleCDSandProtAlignments(refseq, seqs, needlecmd)
    print "Successfully aligned all sequences."
    print "\nWriting the aligned coding DNA sequences to %s." % outfile
    mapmuts.sequtils.WriteFASTA(cds_alignments, outfile)
    print "\nWriting the aligned protein sequences to %s." % outprotfile
    mapmuts.sequtils.WriteFASTA(prot_alignments, outprotfile)


main() # run the program
