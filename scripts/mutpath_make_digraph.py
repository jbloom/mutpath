#!python

"""Script for building mutational trajectory directed graphs.

Using the mutational paths file written by ``mutpath_get_paths.py``, this 
script builds a directed graph through sequence space showing the mutational
trajectory. This path can be visualized using GraphViz.

Written by Jesse Bloom.
"""


import re
import os
import sys
import mutpath.io
import mutpath.sequtils
import mutpath.plot
import mutpath.trajectory


def main():
    """Main body of script."""
    #
    # output is written to out, currently set to standard out
    out = sys.stdout
    statusinterval = 100 # print update after processing this many trees
    credibleinterval = 0.9 # credible interval for dates and persistence
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
    mutpathfile = mutpath.io.ParseFileList(d, 'mutpathfile')
    if len(mutpathfile) != 1:
        raise ValueError("mutpathfile failed to specify exactly one file")
    mutpathfile = mutpathfile[0]
    nodenamefile = mutpath.io.ParseStringValue(d, 'nodenamefile').strip()
    if nodenamefile == 'None':
        nodenames = {}
    else:
        if not os.path.isfile(nodenamefile):
            raise IOError("Cannot find nodenamefile of %s" % nodenamefile)
        nodenames = mutpath.sequtils.Read(nodenamefile)
        nodenames_d = dict([(seq, head.strip()) for (head, seq) in nodenames])
        if len(nodenames) != len(nodenames_d):
            raise ValueError("nodenames file contains duplicate headers")
        nodenames = nodenames_d
        # add stop codon signs to sequence not ending with stop codons to avoid problems due to translateseqs
        for seq in nodenames.keys():
            if seq[-1] != '*':
                nodenames["%s*" % seq] = nodenames[seq]
    dotfile = mutpath.io.ParseStringValue(d, 'dotfile')
    lasttipdate = mutpath.io.ParseFloatValue(d, 'lasttipdate')
    mutationdates = mutpath.io.ParseStringValue(d, 'mutationdates')
    persistencefile = mutpath.io.ParseStringValue(d, 'persistencefile')
    minweight = mutpath.io.ParseFloatValue(d, 'minweight')
    if not 0 <= minweight <= 1:
        raise ValueError("minweight not between zero and one")
    labelcutoff = mutpath.io.ParseFloatValue(d, 'labelcutoff')
    translateseqs = mutpath.io.ParseBoolValue(d, 'translateseqs')
    out.write("\nBuilding mutational trajectory from mutpathfile %s...\n" % mutpathfile)
    out.flush()
    t = mutpath.trajectory.Trajectory(mutpathfile, translateseqs=translateseqs, printprogress=500)
    out.write("Built trajectory from %d individual paths.\n" % t.npaths)
    n_nodesmin = len([node for (node, f) in t.nodes.iteritems() if f >= minweight])
    n_edgesmin = len([edge for (edge, f) in t.edges.iteritems() if f >= minweight])
    out.write("This trajectory contains %d total nodes, of which %d have weights of at least %.5f.\n" % (len(t.nodes), n_nodesmin, minweight))
    out.write("This trajectory contains %d total edges, of which %d have weights of at least %.5f.\n" % (len(t.edges), n_edgesmin, minweight))
    out.write('\nNow writing the dotfile containing the trajectory as a GraphViz directed graph to %s\n' % dotfile)
    mutpath.trajectory.WriteGraphVizTrajectory(t, dotfile, minweight, labelcutoff, nodenames=nodenames)
    mutationdatesfile = "%s.txt" % mutationdates
    out.write("\nWriting the mutation dates in text format to %s" % mutationdatesfile)
    if mutpath.plot.PylabAvailable():
        mutationdatesplot = "%s.pdf" % mutationdates
        out.write("\nPlotting the mutation dates to %s\n" % mutationdatesplot)
    else:
        out.write("\nMutation dates will not be plotted because pylab is not available\n")
        mutationdatesplot = None
    mutpath.trajectory.WriteMutationDates(t, labelcutoff, credibleinterval, mutationdatesfile, mutationdatesplot, lasttipdate)
    out.write("\nWriting persistence times for specified nodes to %s\n" % persistencefile)
    mutpath.trajectory.WriteNodePersistence(t, nodenames, credibleinterval, persistencefile, labelcutoff)
    out.write('\nScript is complete.\n')
    

if __name__ == '__main__':
    main() # run the script
