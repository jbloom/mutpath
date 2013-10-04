#!python

"""Script for extracting mutational paths from BEAST ``.trees`` files.

Extracts mutational paths from BEAST Markov Jumps ``.trees`` files.

Also creates a merged tree without branch annotations.

Written by Jesse Bloom.
"""


import re
import os
import sys
import mutpath.sequtils
import mutpath.io
import mutpath.tree
import mutpath.parse_tree



def main():
    """Main body of script."""
    #
    # output is written to out, currently set to standard out
    out = sys.stdout
    statusinterval = 100 # print update after processing this many trees
    #
    # Read input file and parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the input file name.")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile %s" % infilename)
    d = mutpath.io.ParseInfile(open(infilename))
    out.write('\nRead input argments from %s\n' % infilename)
    out.write('Read the following key / value pairs:\n')
    for (key, value) in d.iteritems():
        out.write("%s %s\n" % (key, value))
    intreefiles = mutpath.io.ParseFileList(d, 'intreefiles')
    burnin = mutpath.io.ParseIntValue(d, 'burnin')
    mergedtreesfile = mutpath.io.ParseStringValue(d, 'mergedtreesfile')
    mutpathsfile = mutpath.io.ParseStringValue(d, 'mutpathsfile')
    fastafile = mutpath.io.ParseFileList(d, 'fastafile')
    if len(fastafile) != 1:
        raise IOError("More than one fastafile specified")
    fastafile = fastafile[0]
    headers_seqs = dict([(header, seq.upper()) for (header, seq) in mutpath.sequtils.Read(fastafile)])
    seqtype = mutpath.io.ParseStringValue(d, 'seqtype').upper()
    if seqtype == 'DNA':
        allowambiguousmismatches = set(['-', 'N', 'R', 'W', 'Y', 'M', 'K', 'S', 'H', 'B', 'V', 'D'])
    elif seqtype == 'PROTEIN':
        allowambiguousmismatches = set(['-', 'X'])
    else:
        raise ValueError("Invalid seqtype of %s" % seqtype)
    startseq = mutpath.io.ParseStringValue(d, 'startseq')
    endseq = mutpath.io.ParseStringValue(d, 'endseq')
    #
    # begin looping over trees
    # bamatch matches branch mutation annotations
    bamatch = re.compile('\[\&(history\_all|h)\=\{\{\d+\,\d+\.\d+,[A-z],[A-z]\}(\,\{\d+\,\d+\.\d+,[A-z],[A-z]\})*\}\]')
    treefilepreface = '' # information before trees lines
    ntreelines = 0
    out.write("\nAll of the non-burnin trees will be written (without branch mutation annotations) to %s.\n" % mergedtreesfile)
    out.write("\nAll of the non-burnin mutational paths will be written to %s.\n" % mutpathsfile)
    trees_out = open(mergedtreesfile, 'w')
    paths_out = open(mutpathsfile, 'w')
    ipath = 0
    for intreefile in intreefiles:
        out.write("\nReading trees from %s...\n" % intreefile)
        numbers_to_strains = mutpath.parse_tree.ReadTreeTranslateBlock(intreefile)
        assert len(numbers_to_strains) == len(headers_seqs), "%d in numbers_to_strains, %d in headers_seqs" % (len(numbers_to_strains), len(headers_seqs))
        tipnames_dict = {}
        for (number, strain) in numbers_to_strains.iteritems():
            try:
                seq = headers_seqs[strain]
            except KeyError:
                raise ValueError("fastafile did not specify sequence for %s" % strain)
            tipnames_dict[number] = {'strain':strain, 'sequence':seq}
        treestarted = fileended = False
        preface = []
        itree = 0
        for line in open(intreefile):
            if fileended:
                raise IOError("file %s has and 'End;' line before the end." % intreefile)
            if (not treestarted) and line[ : 10] == 'tree STATE':
                treestarted = True
            if treestarted:
                if line[ : 4] == 'End;':
                    fileended = True
                    continue
                itree += 1
                if itree < burnin:
                    pass # burn in
                elif itree == burnin:
                    out.write("Finished reading %d burnin trees from %s...\n" % (burnin, intreefile))
                    out.flush()
                else:
                    if (itree - burnin) % statusinterval == 0:
                        out.write("Processed %d trees from %s...\n" % (itree - burnin, intreefile))
                        out.flush()
                    ntreelines += 1
                    trees_out.write(re.sub(bamatch, '', line))
                    # get the Newick string with mutations annotated by history_all key in branch_info 
                    newickstring = mutpath.parse_tree.GetTreeString(line.replace('&h=', '&history_all='))
                    # create tree with 'strain' and 'sequence' keys in info attribute
                    t = mutpath.tree.Tree(newickstring, tipnames_dict=tipnames_dict)
                    mutpath.parse_tree.AssignMutations(t.GetRoot())
                    mutpath.parse_tree.BuildSeqsFromMutsAndTips(t, allowambiguousmismatches)
                    mutpath.tree.AssignTimes(t)
                    path = mutpath.parse_tree.TraceMutationPath(t, startseq, endseq)
                    ipath += 1
                    paths_out.write("MUTPATH %d\n%s\n\n" % (ipath, path))
            else:
                preface.append(line)
                if not treefilepreface:
                    trees_out.write(line)
        preface = ''.join(preface)
        if treefilepreface:
            if treefilepreface != preface:
                raise IOError("The preface lines in %s do not match that from the first *.tree file" % intreefile)
        else:
            treefilepreface = preface
        if not fileended:
            raise IOError("Failed to find an 'End;' line at the end of %s" % intreefile)
    out.write("\nRead a total of %d trees from the %d intreefiles after removing the first %d burnin trees from each file.\n" % (ntreelines, len(intreefiles), burnin))

    trees_out.write('End;')
    trees_out.close()
    paths_out.close()
    out.write("\nWrote %d trees to mergedtreesfile %s." % (ntreelines, mergedtreesfile))
    out.write("\nWrote %d paths to mutpathsfile %s." % (ntreelines, mutpathsfile))
    out.write('\nScript is complete.\n')

    

if __name__ == '__main__':
    main() # run the script
