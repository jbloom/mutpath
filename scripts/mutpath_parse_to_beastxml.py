#!python

"""Script for parsing sequences and building BEAST input XML file.

Written by Jesse Bloom.
"""


import re
import os
import sys
import mutpath.sequtils
import mutpath.io
import mutpath.align
import mutpath.beast
import mutpath.plot


def main():
    """Main body of script."""
    #
    # output is written to out, currently set to standard out
    out = sys.stdout
    #
    # Read input file and parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the input file name.")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile %s" % infilename)
    d = mutpath.io.ParseInfile(open(infilename))
    out.write("\nRead input arguments from %s\n" % infilename)
    out.write("Read the following key / value pairs:\n")
    for (key, value) in d.iteritems():
        out.write("%s %s\n" % (key, value))
    #
    # Now get the sequences
    musclepath = d['musclepath'].strip()
    if not os.path.isfile(d['startseq']):
        raise IOError("Cannot find startseq file %s" % d['startseq'])
    startseq = mutpath.sequtils.Read(d['startseq'])
    if len(startseq) != 1:
        raise IOError("Failed to read exactly one sequence from startseq file %s" % d['startseq'])
    startseq = startseq[0]
    out.write("\nRead the following startseq from %s\n>%s\n%s\n" % (d['startseq'], startseq[0], startseq[1]))
    if not os.path.isfile(d['endseq']):
        raise IOError("Cannot find endseq file %s" % d['endseq'])
    endseq = mutpath.sequtils.Read(d['endseq'])
    if len(endseq) != 1:
        raise IOError("Failed to read exactly one sequence from endseq file %s" % d['endseq'])
    endseq = endseq[0]
    out.write("\nRead the following endseq from %s\n>%s\n%s\n" % (d['endseq'], endseq[0], endseq[1]))
    if not os.path.isfile(d['seqfile']):
        raise IOError("Cannot find seqfile file %s" % d['seqfile'])
    seqs = mutpath.sequtils.Read(d['seqfile'])
    out.write("\nRead %d sequences from seqfile file %s:\n" % (len(seqs), d['seqfile']))
    for (head, seq) in seqs:
        out.write(">%s\n%s\n" % (head, seq))
    # 
    # Convert to upper case, maybe translate, maybe purge redundant
    startseq = (startseq[0], startseq[1].upper())
    endseq = (endseq[0], endseq[1].upper())
    a = mutpath.align.Align([startseq, endseq], musclepath)
    a = mutpath.align.StripGapsToFirstSequence(a)
    a = mutpath.align.StripLeadingTrailingGapsToFirstSequence(a)
    endseq = a[1]
    assert len(startseq[1]) == len(endseq[1])
    for i in range(len(seqs)):
        seqs[i] = (seqs[i][0], seqs[i][1].upper())
    translateseqs = mutpath.io.ParseBoolValue(d, 'translateseqs')
    if translateseqs:
        out.write("\ntranslateseqs is True, so translating all sequences.\n")
        try:
            startseq = mutpath.sequtils.Translate([startseq])[0]
        except:
            sys.stderr.write("\nError translating startseq:\n>%s\n%s\n" % (startseq))
            raise
        try:
            endseq = mutpath.sequtils.Translate([endseq])[0]
        except:
            sys.stderr.write("\nError translating endseq:\n>%s\n%s\n" % (endseq))
            raise
        for i in range(len(seqs)):
            (head, seq) = seqs[i]
            try:
                (head, seq) = mutpath.sequtils.Translate([(head, seq)], truncate_incomplete=True, readthrough_n=True)[0]
            except:
                sys.stderr.write("\nError translating this sequence:\n>%s\n%s\n" % (head, seq))
                raise
            seqs[i] = (head, seq)
    else:
        out.write("\ntranslateseqs is False, so not translating sequences.\n")
    seqs = [startseq, endseq] + seqs # combine all sequences
    purgeredundant = mutpath.io.ParseBoolValue(d, 'purgeredundant')
    if purgeredundant:
        out.write("\npurgeredundant is True, so purging redundant sequences.\n")
        nstart = len(seqs)
        seqs = mutpath.sequtils.PurgeDuplicates(seqs)
        out.write("After purging duplicates, %d of %d sequences remain.\n" % (len(seqs), nstart))
    else:
        out.write("\npurgeredundant is False, so not purging redundant sequences.\n")
    #
    # excludeseqs may specify the manual removal of some sequences
    excludeseqs = d['excludeseqs'].strip()
    if excludeseqs == 'None':
        out.write('\nexcludeseqs set to None, so no sequences specified for exclusion.\n')
    else:
        if not os.path.isfile(excludeseqs):
            raise IOError("Failed to find excludeseqs file of %s" % excludeseqs)
        to_exclude = dict([(line.strip(), True) for line in open(excludeseqs).readlines() if not line.isspace()])
        if not to_exclude:
            out.write('\nexcludeseqs of %s is an empty file, so no sequences specified for exclusion.\n' % excludeseqs)
        else:
            out.write('\nexcludeseqs of %s specifies the following sequences for removal:\n' % excludeseqs)
            retainedseqs = []
            for (head, seq) in seqs:
                for excludehead in to_exclude.keys():
                    if head in excludehead:
                        out.write("excluding: %s\n" % head)
                        break
                else:
                    head = head.replace('[', '(') # square brackets make TreeAnnotator crash
                    head = head.replace(']', ')')
                    retainedseqs.append((head, seq))
            seqs = retainedseqs
            out.write('After removing specified sequences, %d total remain.\n' % len(seqs))
    #
    # align the sequences
    identities = []
    gaps = []
    if not os.path.isdir(musclepath):
        raise IOError("musclepath not found: %s" % musclepath)
    out.write("\nDoing individual pairwise alignments of all seqs to startseq.\n")
    out.write("Gaps relative to startseq are stripped from sequences.\n")
    out.write("Performing alignments with MUSCLE in %s.\n" % musclepath)
    alignedseqs = [seqs[0]]
    assert alignedseqs[0] == startseq
    for iseq in range(1, len(seqs)):
        if iseq % 50 == 0:
            out.write("Completed alignment %d of %d.\n" % (iseq, len(seqs)))
        a = mutpath.align.Align([startseq, seqs[iseq]], musclepath)
        print a
        a = mutpath.align.StripGapsToFirstSequence(a)
        print a
        a = mutpath.align.StripLeadingTrailingGapsToFirstSequence(a)
        print a
        (ident, gap) = mutpath.align.PairwiseStatistics(a)
        identities.append((ident, a[1][0]))
        gaps.append((gap, a[1][0]))
        alignedseqs.append(a[1])
        assert len(a[1][1]) == len(startseq[1]) == len(endseq[1]), "%d, %d, %d, %s" % (len(a[1][1]), len(startseq[1]), len(endseq[1]), a[1][0])
    assert len(alignedseqs) == len(seqs)
    seqs = alignedseqs
    out.write("All %d sequences have been aligned.\n" % len(seqs))
    #
    # make divergence rankings
    identities.sort()
    gaps.sort()
    gaps.reverse()
    divergenceranking = d['divergenceranking'].strip()
    if divergenceranking != 'None':
        gapsfile = "%s_gaps.txt" % divergenceranking
        identitiesfile = "%s_identities.txt" % divergenceranking
        out.write("\nWriting information about gaps and identities (for identifying potential outliers) to %s and %s.\n" % (gapsfile, identitiesfile))
        for (fname, xlist, header) in [(gapsfile, gaps, 'fraction_gaps'), (identitiesfile, identities, 'fraction_identities')]:
            f = open(fname, 'w')
            f.write("#name\t%s\n" % header)
            for (x, head) in xlist:
                f.write("%s\t%.5f\n" % (head, x))
            f.close()
        if mutpath.plot.PylabAvailable():
            gapsplot = "%s_gaps.pdf" % divergenceranking
            identitiesplot = "%s_identities.pdf" % divergenceranking
            out.write("Using pylab to make cumulative distribution plots %s and %s" % (gapsplot, identitiesplot))
            mutpath.plot.CumulativeFractionPlot([x[0] for x in gaps], gapsplot, 'gaps per sequence', 'fraction of sites that are gaps')
            mutpath.plot.CumulativeFractionPlot([x[0] for x in identities], identitiesplot, 'fraction identities among aligned positions', 'fraction of identities')
        else:
            out.write("Not making any plots as pylab is not available.\n")
    #
    # write output files
    fastaoutfile = "%s.fasta" % d['outfileprefix']
    out.write("\nNow writing these sequences in FASTA format to %s\n" % fastaoutfile)
    headmatch = re.compile('^\S+ +(?P<strain>[\S ]+) +(?P<date>\d+/\d*/\d*|\d+\-\d*\-\d*) +\S+$')
    n_with_name = {}
    for i in range(len(seqs)):
        (head, seq) = seqs[i]
        m = headmatch.search(head)
        if not m:
            raise ValueError("Failed to match header:\n%s" % head)
        strain = m.group('strain').strip().replace(' ', '_')
        date = mutpath.sequtils.DateToOrdinal(m.group('date'), 1968)
        newhead = "%s_%.2f" % (strain, date)
        j = 1 # keep track of multiple sequences with the same name
        while newhead in n_with_name:
            j += 1
            newhead = '%s_%d_%.2f' % (strain, j, date)
        n_with_name[newhead] = True
        seqs[i] = (newhead, seq)
    mutpath.sequtils.Write(seqs, fastaoutfile)
    nexusoutfile = '%s.nex' % d['outfileprefix']
    out.write('\nNow writing these sequences in NEXUS format to %s\n' % nexusoutfile)
    if translateseqs:
        seqtype = 'PROTEIN'
    else:
        # autodetect sequence type
        seqtype = 'DNA'
        for nt_or_aa in seqs[0][1]:
            if nt_or_aa.upper() not in ['A', 'T', 'C', 'G', 'N']:
                seqtype = 'PROTEIN'
                break
    mutpath.sequtils.WriteNEXUS(seqs, nexusoutfile, seqtype)
    sitemodel = d['sitemodel']
    if sitemodel == 'HKY' and translateseqs:
        raise ValueError('set sitemodel to HKY (nucleotides) but translateseqs to True')
    if sitemodel == 'HKY':
        datatype = 'nucleotide'
    elif sitemodel == 'JTT':
        datatype = 'amino acid'
    else:
        raise ValueError("Invalid sitemodel of %s" % sitemodel)
    xmloutfile = "%s.xml" % d['outfileprefix']
    out.write('\nNow writing these sequences into a BEAST xml file of %s\n' % xmloutfile)
    usemarkovjumps = mutpath.io.ParseBoolValue(d, 'usemarkovjumps')
    chainlength = mutpath.io.ParseIntValue(d, 'chainlength')
    screenloginterval = mutpath.io.ParseIntValue(d, 'screenloginterval')
    fileloginterval = mutpath.io.ParseIntValue(d, 'fileloginterval')
    mutpath.beast.WriteXML(xmloutfile, seqs, datatype, chainlength, screenloginterval, fileloginterval, usemarkovjumps, sitemodel)
    out.write('\nScript is complete.\n')

    

if __name__ == '__main__':
    main() # run the script
