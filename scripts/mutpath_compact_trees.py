#!python

"""Script for compacting BEAST trees run with Markov Jumps mutation mapping.

Written by Jesse Bloom.
"""


import re
import os
import sys


def RoundString(m, decimaldigits):
    """
    Takes a string for a float and returns version with specified decimal digits.

    *m* is an *re* match object. The first group (*m.group(0)*)
    must represent a string giving a number that can be
    converted to a float using *float*.

    *decimaldigits* is an integer >= 0 that represents how
    many digits we retain in the rounded number.

    Example:

    >>> num_match = re.compile('\d+\.\d+')
    >>> m = num_match.search('53.2345145')
    >>> RoundString(m, 4)
    '53.2345'

    """
    if not (isinstance(decimaldigits, int) and decimaldigits >= 0):
        raise ValueError("decimaldigits not an integer >= 0: %s" % str(decimaldigits))
    x = float(m.group(0))
    mstring = '%.' + str(decimaldigits) + 'f'
    return mstring % x



def main():
    """Main body of script."""
    #
    # number of digits we retain in decimals
    decimaldigits = 4
    def RoundToDD(m):
        """Uses RoundString with decimaldigits digits"""
        return RoundString(m, decimaldigits)
    #
    # output is written to out, currently set to standard out
    out = sys.stdout
    #
    # Read input file and parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the input *.trees file name.")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile %s" % infilename)
    (base, ext) = os.path.splitext(infilename)
    if ext != '.trees':
        raise ValueError('Input tree file of %s does not have the extension *.trees.' % infilename)
    outfilename = "%s_compact%s" % (base, ext)
    out.write("\nWe will compact the file %s to make the more compact version %s.\n" % (infilename, outfilename))
    if os.path.isfile(outfilename):
        raise IOError("Output file of %s already exists. We will not overwrite it." % outfilename)
    out.write('\nNow compacting the file...')
    out.flush()
    replace_with_null = [ # patterns to remove from node and branch annotations
        re.compile('c\_allTransitions\[\d+\]\=\d+\.\d+\,'),
        re.compile('c\_allTransitions\[\d+\]\=\d+\.\d+'),
        re.compile('\&states="[A-Z]+"\,'),
        re.compile('\&states="[A-Z]+"'),
        re.compile('states="[A-Z]+"'),
        re.compile('\[\]'),
        ]
    history_match = re.compile('history\_all\=\{')
    float_match = re.compile('\d+\.\d{%d,}' % decimaldigits)
    try:
        fin = open(infilename)
        fout = open(outfilename, 'w')
        iline = 0
        for line in fin:
            if line[ : 5] != 'tree ':
                fout.write(line) # not a trees line, so do not change
                continue
            for m in replace_with_null:
                line = re.sub(m, '', line)
            line = re.sub(history_match, 'h={', line)
            line = re.sub(float_match, RoundToDD, line)
            fout.write(line)
            iline += 1
            if iline % 50 == 0:
                out.write("\nCompacted %d tree lines..." % iline)
                out.flush()
    finally:
        fin.close()
        fout.close()
    out.write('\nCompleted creating the compacted file %s.\n' % outfilename)
    out.write('\nThe original size of %s is %d bytes; the compacted version in %s is %d bytes.\n' % (infilename, os.path.getsize(infilename), outfilename, os.path.getsize(outfilename)))
    out.write('\nYou are now free to manually delete the larger %s file.\n' % infilename)
    out.write('\nScript is complete.\n')

    

if __name__ == '__main__':
    main() # run the script
