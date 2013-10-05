#!python

"""Script for annotating a *.trees file to be read by FigTree.

Written by Jesse Bloom.
"""


import re
import os
import sys
import mutpath.io



def main():
    """Main body of script."""
    #
    # number of digits we retain in decimals
    out = sys.stdout
    #
    # Read input file and parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the input file name.")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile %s" % infilename)
    lines = [line for line in open(infilename).readlines() if line[0] != '#' and not line.isspace()]
    intreefile = outtreefile = None
    name_d = {}
    for line in lines:
        entries = line.split(None, 1)
        if len(entries) != 2:
            raise ValueError("line does not contain two entries:\n" % line)
        if entries[0].strip() == 'intreefile':
            if intreefile == None:
                intreefile = entries[1].strip()
                out.write("\nThe intreefile is %s" % intreefile)
                if not os.path.isfile(intreefile):
                    raise IOError("Cannot find intreefile of %s" % intreefile)
            else:
                raise ValueError("duplicate entries for intreefile")
        elif entries[0].strip() == 'outtreefile':
            if outtreefile == None:
                outtreefile = entries[1].strip()
                out.write("\nThe outtreefile is %s" % intreefile)
            else:
                raise ValueError("duplicate entries for outtreefile")
        else:
            (name, abbr) = (entries[0].strip(), entries[1].strip())
            if name in name_d:
                raise ValueError("duplicate entries for tip name %s" % name)
            out.write("\nTip %s will be renamed to %s" % (name, abbr))
            name_d[name] = abbr
    startedtaxa = finishedtaxa = False
    try:
        fin = open(intreefile)
        fout = open(outtreefile, 'w')
        for line in fin:
            newline = line
            line = line.strip()
            if startedtaxa and not finishedtaxa:
                if line == ';':
                    finishedtaxa = True
                else:
                    if line[1 : -1] in name_d:
                        newline = '%s[&!name="%s"]\n' % (newline[ : -1], name_d[line[1 : -1]])
                        del name_d[line[1 : -1]]
                    else:
                        newline = '%s[&!name=""]\n' % (newline[ : -1])
            elif (line.lower() == 'taxlabels'):
                if startedtaxa or finishedtaxa:
                    raise IOError("Found Taxlabels line but already started or finished taxa block")
                startedtaxa = True
            else:
                pass
            fout.write(newline)
        if not (startedtaxa and finishedtaxa):
            raise IOError("Failed to find start and end of Taxlabels taxa block")
        if name_d:
            raise IOError("Failed to find the following taxa:\n%s" % '\n'.join(name_d.values()))
        fout.write('\n\nbegin figtree;')
        fout.write('\n  set tipLabels.displayAttribute="Names";')
        fout.write('\n  set trees.order=true;')
        fout.write('\n  set tipLabels.fontSize=12;')
        fout.write('\n  set trees.orderType="decreasing";')
        fout.write('\nend;\n')
    finally:
        fin.close()
        fout.close()
    out.write('\nCompleted creating outtreefile %s.\n' % outtreefile)
    out.write('\nScript is complete.\n')

    

if __name__ == '__main__':
    main() # run the script
