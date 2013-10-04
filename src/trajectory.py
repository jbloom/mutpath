"""Module for representing mutational trajectories as directed graphs.

Represents mutational trajectories through sequence space, which is the space in
which each node is a unique sequence and edges are directional connections
between nodes corresponding to mutations.

These digraphs can be used to visualize mutational trajectories through sequence
space. They are designed to be visualized using the GraphViz program.

Written by Jesse Bloom.


Functions defined in this module
------------------------------------
`WriteGraphVizTrajectory` - Writes a GraphViz visualization of a *Trajectory* object.

`WriteMutationDates` - Writes files giving mutations dates and credible intervals.

`WriteNodePersistence` - Writes times that nodes persist (time before next mutation).

`DistanceAlongPath` - Distance along a mutational path.

`HeuristicTraceback` - Tries to find the last high weight predecessor of a node.

`IteratePaths` - Iterates of paths in a mutational path file.


Classes defined in this module
--------------------------------
`Trajectory` - A class for representing a mutational trajectory through sequence space.


Detailed documentation for functions
--------------------------------------
Provided in the individual function documentation strings below.
"""


import os
import re
import sys
import math
import sequtils
import stats
import plot


def DistanceAlongPath(startseq, endseq, s):
    """Returns the distance along the mutational path.

    This distance of a sequence *s* along the path is defined 
    as the Hamming Distance between *startseq* and *s* minus
    the Hamming Distance between *endseq* and *s* plus the
    Hamming distance between *startseq* and *endseq*.

    Sequences are not case sensitive.
    """
    assert len(s) == len(startseq) == len(endseq)
    s_to_start = len([1 for (x, y) in zip(s.upper(), startseq.upper()) if x != y])
    s_to_end = len([1 for (x, y) in zip(s.upper(), endseq.upper()) if x != y])
    start_to_end = len([1 for (x, y) in zip(startseq.upper(), endseq.upper()) if x != y])
    return s_to_start - s_to_end + start_to_end



def HeuristicTraceback(t, node, cutoff):
    """Traces back to find last high weight precessor of a node.

    *t* is a *Trajectory* object.

    *node* is the string name of a node in *t*.

    *cutoff* is a number > 0 and <= 1.

    This function starts at node, and traces back along the trajectory
    to return the first predecessor of *node* with a weight >= *cutoff*.
    It does this by tracing back from *node* along its highest weight
    incoming edge to that predecessor, and then from that predecessor
    along its highest weight edge, etc until we find a predecessor
    with weight >= *cutoff*. This approach
    is not absolutely guaranteed to find the first predecessor with weight
    > *cutoff*, but it should work for reasonable trajectories. But beware,
    hence the word 'Heuristic' in the name of this function.

    The return value is the string name of the first high weight predecessor.

    This function recursively calls itself.
    """
    assert node in t.nodes
    weights_predecessors = []
    for ((n1, n2), weight) in t.edges.iteritems():
        if n2 == node:
            weights_predecessors.append((weight, n1))
    if not weights_predecessors:
        raise ValueError("failed to find predecessor")
    weights_predecessors.sort()
    weights_predecessors.reverse()
    (weight, predecessor) = weights_predecessors[0]
    if t.nodes[predecessor] >= cutoff:
        return predecessor
    else:
        return HeuristicTraceback(t, predecessor, cutoff)


def WriteNodePersistence(t, nodestowrite, interval, persistencefile, cutoff):
    """Writes times for which nodes persist before the next mutation.

    The trajectory *t* specifies a set of nodes. For each node specified
    by *nodestowrite* and with weight >= *cutoff*, reports the time until that 
    node experiences the next mutation that moves it to a new node.
    If *t* is a trajectory through protein sequence space that was
    created from nucleotide sequences (this will be the case if *t*
    was created with *translateseqs = True*), the next mutation that
    moves it to a new node is a nonsynonymous mutation. This method
    then also records the time until the first mutation of any type
    (synonymous or nonsynonymous) after the trajectory moves to nodes.
    The persistence is written as the posterior median from all 
    paths containing plus the Bayesian credible interval specified 
    by *interval*. The persistence
    times are written to the text file *persistencefile*.

    CALLING VARIABLES:

    * *t* is a *Trajectory* object that contains the persistence data.
    
    * *nodestowrite* is a dictionary specifying for which nodes we
      write persistence times, and the names used when writing these
      nodes. It is keyed by node sequences (which
      are the identifiers for the nodes in t.nodes) and the values are
      strings giving names that are used to label the nodes in the 
      output. However, there does not actually have to be a
      node with persistence data for each key in *nodestowrite* -- if there
      is not persistence data for a node key, nothing is written for it.

    * *interval* specifies the range of the Bayesian credible interval,
      for example a value of 0.9 means that we print the 90% credible
      interval.

    * *persistencefile* is a string giving the name of the text file
      that we create which contains the persistence times. It has
      headers explaning its format. If this file already exists, it
      is overwritten. The order in which nodes is written is arbitrary.

    * *cutoff* is a weight cutoff (fraction of paths containing this node).
      We only write persistence times for nodes that are both in 
      *nodestowrite* and have weights >= *cutoff*. This keeps us from
      writing persistence times for nodes that only occur rarely.

    """
    d = dict([(name, node) for (node, name) in nodestowrite.iteritems()])
    f = open(persistencefile, 'w')
    f.write("# Name : name of the node\n")
    f.write("# MedianPersistence : posterior median time to next node\n")
    f.write("# MinPersistenceInterval : minimum of %.2f percent credible interval\n" % (interval * 100))
    f.write("# MaxPersistenceInterval : maximum of %.2f percent credible interval\n" % (interval * 100))
    if t.persistence != t.timetofirstmut:
        f.write("# MedianTimeToFirstMut : posterior median time to first mutation\n")
        f.write("# MinTimeToFirstMutInterval : minimum of %.2f percent credible interval\n" % (interval * 100))
        f.write("# MaxTimeToFirstMutInterval : maximum of %.2f percent credible interval\n" % (interval * 100))
        f.write("#\n")
        f.write("#Name\tMedianPersistence\tMinPersistenceInterval\tMaxPersistenceInterval\tMedianTimeToFirstMut\tMinTimeToFirstMut\tMaxTimeToFirstMut\n")
    else:
        f.write("#\n")
        f.write("#Name\tMedianPersistence\tMinPersistenceInterval\tMaxPersistenceInterval\n")
    for (node, persistence) in t.persistence.iteritems():
        if (node not in nodestowrite):
            continue
        if (len(persistence) / float(t.npaths)) < cutoff:
            continue
        name = nodestowrite[node]
        (median, mininterval, maxinterval) = stats.MedianCredibleInterval(persistence, interval)
        f.write("%s\t%f\t%f\t%f" % (name, median, mininterval, maxinterval))
        if t.persistence != t.timetofirstmut:
            (median, mininterval, maxinterval) = stats.MedianCredibleInterval(t.timetofirstmut[node], interval)
            f.write("\t%f\t%f\t%f\n" % (median, mininterval, maxinterval))
        else:
            f.write("\n")
    f.close()


def WriteMutationDates(t, labelcutoff, interval, datesfile, datesplot, lasttipdate):
    """Creates text file and plot showing dates of mutations.

    For each mutation that occurs in at least *labelcutoff* fraction of 
    the paths that form the trajectory *t*, this function writes the
    posterior median and a Bayesian credible interval for the date of
    first occurrence of that mutation. The posterior is taken over all
    paths that contain that mutation. 

    The output also provides information about whether the mutations
    are on the branch from the starting sequence to the common ancestor,
    or from the common ancestor to the starting sequence.

    CALLING VARIABLES:

    * *t* is a *Trajectory* object that contains the mutations data.

    * *labelcutoff* is a number > 0 and <= 1. We write the dates of all
      mutations in *t* that have weights >= *labelcutoff* (occur in
      at least this fraction of the paths).

   * *interval* specifies the range of the Bayesian credible interval,
     for example a value of 0.9 means that we print the 90% credible
     interval.

    * *datesfile* is the name of the text file that we create which
      contains the dates and intervals. It is overwritten if it does
      not already exist.

    * *datesplot* is the name of the plot file that we create using
      matplotlib. This plot can only be created if matplotlib is 
      available. So first check on this (for example using 
      *Plot.PylabAvailable()*. If matplotlib is not available, or if
      you don't want to make the plot, make this argument *None*. Otherwise
      make it the name of the PDF plot file that you want to create.

    * *lasttipdate* specifies the absolute units for the dates. Dates
      in *t* will be specified in units of time before the most recent
      tip. Here provide a number giving the date of the most recent tip,
      and the dates shown are then this number minus the time for
      each mutation.
    """
    mutdatestoca = [] # keys are (median, mininterval, maxinterval, mut, fractoca, weight)
    mutdatesfromca = [] # keys are (median, mininterval, maxinterval, mut, fractoca, weight)
    n = t.npaths
    for (mut, muttimes) in t.mutations.iteritems():
        nmut = len(muttimes)
        weight = nmut / float(n)
        if weight >= labelcutoff:
            # mutation meets the cutoff
            fractoca = t.mutationstoca[mut] / float(nmut)
            (median, mininterval, maxinterval) = stats.MedianCredibleInterval(muttimes, interval)
            # we interchange minimum and median on the next line because the times
            # are in units before last tip prior to this conversion
            (median, maxinterval, mininterval) = (lasttipdate - median, lasttipdate - mininterval, lasttipdate - maxinterval)
            if fractoca > 0.5:
                mutdatestoca.append((median, mininterval, maxinterval, mut, fractoca, weight))
            else:
                mutdatesfromca.append((median, mininterval, maxinterval, mut, fractoca, weight))
    mutdatestoca.sort()
    mutdatestoca.reverse()
    mutdatesfromca.sort()
    mutdates = mutdatestoca + mutdatesfromca
    f = open(datesfile, 'w')
    f.write('# Mutation : mutation in 1, 2, ... numbering\n')
    f.write('# FracOccurrence : fraction of paths containing this mutation\n')
    f.write('# FracToCommonAncestor : fraction of times which this mutation is on path from starting sequence to common ancestor\n')
    f.write('# MedianDate : posterior median date of mutation\n')
    f.write('# MinInterval : minimum of %.2f percent Bayesian credible interval (median centered)\n' % interval)
    f.write('# MaxInterval : maximum of %.2f percent Bayesian credible interval (median centered)\n' % interval)
    f.write('#\n')
    f.write('# Mutation\tFracOccurrence\tFracToCommonAncestor\tMedianDate\tMinInterval\tMaxInterval\n')
    for (median, mininterval, maxinterval, mut, fractoca, weight) in mutdates:
        f.write('%s\t%f\t%f\t%f\t%f\t%f\n' % (mut, weight, fractoca, median, mininterval, maxinterval))
    f.close()
    if datesplot:
        plot.DatesPlot(mutdates, datesplot, interval)



def WriteGraphVizTrajectory(t, graphvizfile, minweight, labelcutoff,\
        nodenames=None, nodesize=0.9, ranksep=0.1, nodesep=0.2, penwidth=10,\
        fontsize=40, arrowsize=1.4, fontname='Helvetica-Bold', rankdir='TB',\
        startendasdiamonds=True):
    """Writes a GraphViz visualization of a *Trajectory* object.

    This function creates a file *graphvizfile* that can be used to visualize the
    directed graph represented by a *Trajectory* object *t*. Graphviz is a freely
    available software package (http://www.graphviz.org/) for visualizing graphs.
    The trajectory is written in the DOT language 
    (http://www.graphviz.org/doc/info/lang.html).

    The areas of nodes and the widths of edges are proportional to their weights.

    The color saturations of nodes and edges are linearly proportional to their 
    weights.

    The rank of nodes (for example, their vertical position when *rankdir* 
    is 'LR') is ordered according to their distance along the path from the
    starting to ending sequence. This distance is defined as the Hamming
    Distance from the starting node minus the Hamming Distance from the ending node
    plus the Hamming distance between the starting and ending nodes.

    CALLING VARIABLES:

    * *t* is the *Trajectory* object that contains the trajectory that we want
      to visualize.

    * *graphvizfile* is a string giving the name of the GraphViz input file that
      we want to create. It will typically end with the extension ``.dot``. If this
      file already exists, it is overwritten. You should be able to open this file
      directly with Graphviz. The file is written in the DOT language
      (http://www.graphviz.org/doc/info/lang.html).

    * *minweight* is a number specifying the minimum weight that a node or edge
      must possess in order to be shown on the graph. Nodes or edges with 
      weights < *minweight* are not included. Note that this creates a possibility
      for orphan nodes if a node has a weight >= *minweight* but all of its
      incoming and outgoing nodes have weights < *minweight*. To show all nodes
      and edges regardless of weight, set *minweight* to zero. However, this can
      sometimes lead to a very large *graphvizfile* since there can be a huge
      number of very low weight nodes / edges.

    * *labelcutoff* is the minimum weight that an edge must possess in order
      to be labeled on the graph. In addition, all nodes with weight >=
      *labelcutoff* have an incoming edge that is labeled. If there is not
      such an edge, then traces back to find the first predecessor node with
      weight *labelcutoff* and then draws a different colored edge spanning
      multiple mutations to connect these nodes. Generally, you would want
      *labelcutoff > 0.5*.

    * *nodenames* is an optional argument that allows you to specify names
      for nodes. It is *None* by default. If you set it to another value,
      it should be a dictionary. It is keyed by node sequences (which
      are the identifiers for the nodes in t.nodes) and the values are
      strings giving names that are used to label the nodes in the 
      trajectory. These names are written on the nodes only if
      the weight for that node is >= *labelcutoff*.

    OPTIONAL CALLING VARIABLES SPECIFYING FORMATTING DETAILS:

    * *nodesize* is the height of a node with weight.

    * *ranksep* is the separation between ranks, as fraction of *nodesize*.

    * *nodesep* is the minimum separation between nodes of the same rank,
      as fraction of *nodesize*.

    * *penwidth* is the pen width of an edge.

    * *fontsize* is the font size.

    * *arrowsize* is the size of the arrows.

    * *fontname* is the font style.

    * *rankdir* specifies the direction the ranks move. If set to 'TB' then 
      the graph moves from top to bottom. If set to 'LR' then the graph moves
      from left to right.

    * *startendasdiamonds* is a Boolean switch. If True, we show the starting
      and ending nodes as diamonds rather than circles. We also make these
      starting and ending nodes larger in size to fit their full labels. If 
      False, we make them circles with size proportional to weights like
      all other nodes.

    """
    f = open(graphvizfile, 'w')
    f.write('digraph G { rankdir=%s; ranksep=%f; nodesep=%f;\n' % (rankdir, ranksep * nodesize, nodesep * nodesize))
    # first write the nodes ordered into subgraphs of the same rank by DistanceAlongPath
    nodes_by_d = {}
    needs_incoming = {} # does node need an incoming edge?
    for (node, weight) in t.nodes.iteritems():
        if weight < minweight:
            continue # weight too low
        d = DistanceAlongPath(t.startseq, t.endseq, node)
        if d in nodes_by_d:
            nodes_by_d[d].append((node, weight))
        else:
            nodes_by_d[d] = [(node, weight)]
        if (weight >= labelcutoff) and node != t.startseq:
            needs_incoming[node] = True
    for d in range(max(nodes_by_d.keys()) + 1):
        if d not in nodes_by_d:
            continue # none of this distance
        f.write('\tsubgraph %d { label="DistanceAlongPath%d"; rank=same;\n' % (d, d))
        for (node, weight) in nodes_by_d[d]:
            if startendasdiamonds and (node == t.startseq or node == t.endseq):
                shape = 'diamond'
                fixedsize = 'false'
            else:
                shape = 'circle'
                fixedsize = 'true'
            if nodenames and (node in nodenames) and weight >= labelcutoff:
                nodelabel = "%s" % nodenames[node]
            else:
                nodelabel = ''
            f.write('\t\tnode [style=filled shape=%s label="%s" height=%f color="0.7 %f 0.9" penwidth=%f arrowsize=%f fontsize=%d fontname="%s" fontcolor="white" fixedsize=%s] "%s";\n' % (shape, nodelabel, nodesize * math.sqrt(weight), weight, penwidth, arrowsize, fontsize, fontname, fixedsize, node))
        f.write('\t}\n')
    # now write all of the edges
    # In order to get good overlay, first we write unabeled edges, then
    # labeled edges, and finally implied connections between major nodes without
    # connecting labeled edges.
    labeled_edges = []
    for ((node1, node2), weight) in t.edges.iteritems():
        if weight < minweight:
            continue # weight too low
        if weight >= labelcutoff:
            assert len(node1) == len(node2)
            diffs = [i for i in range(len(node1)) if node1[i] != node2[i]]
            if len(diffs) != 1:
                raise ValueError("Should be exactly one difference")
            i = diffs[0]
            edgelabel = '%s%d%s' % (node1[i], i + 1, node2[i])
            if node2 in needs_incoming:
                del needs_incoming[node2]
        else:
            edgelabel = ''
        edgestring = '\t"%s" -> "%s" [weight=%f penwidth=%f color="0.7 %f 0.9" arrowsize=%f label="%s" fontsize=%d fontname="%s"];\n' % (node1, node2, weight, penwidth * weight, weight, arrowsize, edgelabel, fontsize, fontname)
        if edgelabel:
            labeled_edges.append(edgestring) # write these later
        else:
            f.write(edgestring)
    f.write(''.join(labeled_edges)) # now write labeled edges
    # now find implied connections between major nodes without incoming labeled edges
    for node in needs_incoming:
        predecessor = HeuristicTraceback(t, node, labelcutoff)
        diffs = [i for i in range(len(node)) if node[i] != predecessor[i]]
        assert len(diffs) >= 1
        diffs.sort()
        edgelabel = '-'.join(["%s%d%s" % (predecessor[i], i + 1, node[i]) for i in diffs])
        f.write('\t"%s" -> "%s" [weight=0 penwidth=%f color="0.0 1.0 0.9" arrowsize=%f label="%s" fontsize=%d fontname="%s" fontcolor="0.0 1.0 0.9"];\n' % (predecessor, node, penwidth, arrowsize, edgelabel, fontsize, fontname))
    f.write('}') 
    f.close()



def IteratePaths(pathfile):
    """Iterates over paths in a mutational path file.

    *pathfile* should be a string giving a name of an input file specifying one or
    more mutational paths. These files are of the format created by 
    ``mutpath_get_paths.py``. The required format is detailed below.
   
    This function will iterate over all paths in *pathfile*. For each path,
    it will return the tuple 
    *(startseq, starttime, endseq, endtime, caseq, catime, tocamuts, fromcamuts)*.
    The entries of these tuples are as follows. All sequences are converted
    to upper case, as are all letters in the mutation notations. The times
    are measured in units before the most recent tip of the tree. Tuple entries:

        * *startseq* is the starting sequence specified by *startstrain_seq*

        * *starttime* is the time of *startseq* specified by *startstrain_time*
        
        * *endseq* is the ending sequence specified by *endstrain_seq*

        * *endtime* is the time of *endseq* specified by *endstrain_time*

        * *caseq* is the common ancestor sequence specified by *commonancestor_seq*

        * *catime* is the time of *caseq* specified by *commonancestor_time*

        * *tocamuts* is a list of the mutations going from *startseq* to *caseq*,
          specified in the order they are listed in the file (should be along
          the path) as the 2-tuples of the form *('A6G', 42.713)* where the
          entries are the mutation and then the time.

        * *fromcamuts* is like *tocamuts*, but for mutations going from 
          *caseq* to *endseq*.
    
    The format of *pathfile* is as follows. This file should list
    mutational paths as::

            MUTPATH 1
            startstrain_name A/Aichi/2/1968_1968.50
            startstrain_seq ATGGCAATGGGCTAA 
            startstrain_time 42.5
            endstrain_name A/Brisbane/10/2007_2007.10
            endstrain_seq ATGACGATTGGATAA
            endstrain_time 3.9
            commonancestor_seq ATGGCGATGGGCTAA
            commonancestor_time 43.12713
            startstrain_to_commonancestor_path
            A6G : 42.713     
            commonancestor_to_endstrain_path
            G9T : 31.732
            G4A : 25.1343
            C12A : 10.134

            MUTPATH 2
            startstrain_name A/Aichi/2/1968_1968.50
            startstrain_seq ATGGCAATGGGCTAA 
            startstrain_time 42.5
            endstrain_name A/Brisbane/10/2007_2007.10
            endstrain_seq ATGACGATTGGATAA
            endstrain_time 3.9
            commonancestor_seq ATGGCGATGGGCTAA
            commonancestor_time 44.12713
            startstrain_to_commonancestor_path
            A6G : 42.113     
            G9T : 43.124
            commonancestor_to_endstrain_path
            G4A : 21.1343
            C5A : 19.531
            A5C : 19.402
            C12A : 9.134

    The file lists each of the paths numbered starting at 1. 
    Within each path, the mutations are indicated with numbering starting 
    at 1 for the first position in the sequence. The times for the mutations, 
    the starting and ending strains, and the most recent common ancestor of these 
    two strains, are also indicated. These times are measured in units before 
    the most recent tip node (so the root node would have the largest value of time).
    The mutations must move from the starting to the ending sequence, and if
    multiple paths are specified, then they all must have the same starting and
    ending sequences.
    """
    mutmatch = re.compile('^(?P<mut>[A-z\*\-]\d+[A-z\*\-]) : (?P<time>\d+\.*\d*)$')
    if not os.path.isfile(pathfile):
        raise IOError("Cannot find pathfile %s" % pathfile)
    f = open(pathfile)
    firststartseq = firstendseq = None
    while True:
        try:
            line = f.next()
        except StopIteration:
            break # no more lines
        lines = []
        while not line.isspace():
            lines.append(line.strip())
            line = f.next()
        tocamuts = []
        fromcamuts = []
        assert lines[0][ : 7] == 'MUTPATH'
        assert lines[1][ : 16] == 'startstrain_name'
        assert lines[2][ : 15] == 'startstrain_seq'
        startseq = lines[2].split()[1].strip().upper()
        if firststartseq == None:
            firststartseq = startseq
        elif firststartseq != startseq:
            raise IOError("Change in startseq")
        assert lines[3][ : 16] == 'startstrain_time'
        starttime = float(lines[3].split()[1])
        assert lines[4][ : 14] == 'endstrain_name'
        assert lines[5][ : 13] == 'endstrain_seq'
        endseq = lines[5].split()[1].strip().upper()
        if firstendseq == None:
            firstendseq = endseq
        elif firstendseq != endseq:
            raise IOError("Change in endseq")
        assert lines[6][ : 14] == 'endstrain_time'
        endtime = float(lines[6].split()[1])
        assert lines[7][ : 18] == 'commonancestor_seq'
        caseq = lines[7].split()[1].strip().upper()
        assert lines[8][ : 19] == 'commonancestor_time'
        catime = float(lines[8].split()[1])
        assert lines[9] == 'startstrain_to_commonancestor_path'
        i = 10
        while lines[i] != 'commonancestor_to_endstrain_path' and i < len(lines):
            m = mutmatch.search(lines[i])
            if not m:
                raise ValueError("Failed to match mutation line:\n%s" % lines[i])
            tocamuts.append((m.group('mut'), float(m.group('time'))))
            i += 1
        if i < len(lines):
            if lines[i] != 'commonancestor_to_endstrain_path':
                raise ValueError("Expected 'commonancestor_to_endstrain_path', but got:\n%s" % lines[i])
            i += 1
        while i < len(lines):
            m = mutmatch.search(lines[i])
            if not m:
                raise ValueError("Failed to match mutation line:\n%s" % lines[i])
            fromcamuts.append((m.group('mut'), float(m.group('time'))))
            i += 1
        yield (startseq, starttime, endseq, endtime, caseq, catime, tocamuts, fromcamuts)
    f.close()



class Trajectory(object):
    """Class for representing a mutational trajectory through sequence space.

    This class represents a mutational trajectory in sequence space. The trajectory
    is a directed graph consisting of nodes (sequences) and edges (mutations
    connecting nodes). The trajectory moves from one known sequence to another
    known sequence, passing through some number of uncertain intermediates (nodes).
    The trajectory is created by passing it a set of possible mutational paths
    from the starting to ending sequence. In the trajectory, the weight of each node
    corresponds to the fraction of paths that contain that sequence, while the
    weight of each edge corresponds to the fraction of paths that contain that edge.
    Note that if a path contains a node or edge more than once (which can happen
    if there are mutational cycles), the node or edge is still considered to have
    occurred once in that path for the purposes of assigning the weights.

    Each *Trajectory* object *t* has the following attributes:

        * *t.npaths* : the number of individual paths used to construct the
          overall trajectory.

        * *t.startseq* : a string giving the starting sequence for the trajectory.

        * *t.endseq* : a string giving the ending sequence for the trajectory.

        * *t.nodes* : a dictionary keyed by strings representing the sequences for
          each node found at least once in the trajectory, and with values equal to
          the weight of that node (fraction of paths containing the node).

        * *t.edges* : a dictionary keyed by 2-tuples of strings *(s1, s2)* and
          values giving the weight of the directed edges from sequence *s1* to
          *s2*.

        * *t.mutations* : a dictionary keyed by mutation strings of the
          form 'G5A' specifying mutations where the numbering is in
          1, 2, ... For each mutation that occurs in at least one of
          the paths passed to this trajectory, there will be a key.
          The values are lists giving the times of occurrence for all
          occurrences of that mutation in the paths used to create this
          trajectory. If a mutation occurs more than once in a path,
          only the time for its first occurrence is listed. So the
          fraction of paths that contain some mutation *m* is
          *t.mutations[m] / float(t.npaths)*. Note that if *translateseqs*
          is *True*, then the mutations specified here are only the
          protein mutations, not the nucleotide ones in the underlying
          nucleotide sequences.

        * *t.mutationstoca* : a dictionary keyed by mutation strings just
          as for *t.mutations*. Each mutation that is added to the lists
          in *t.mutations* can arise on the branch from the starting sequence
          to the common ancestor, or on the branch from the common ancestor
          to the ending sequence. The value of *t.mutationstoca[mut]* is
          the number of times that the first occurrence of *mut* is on
          the route from starting sequence to the common ancestor. So if
          *mut* is always on the path from the starting sequence to the
          common ancestor, then *t.mutationstoca[mut] == len(t.mutations[mut])*.
          If it is always on the path from the starting sequence to the 
          ending sequence, then *t.mutations[mut] == 0*.

        * *t.persistence* : a dictionary keyed by the node sequences and
          with the values being a list of numbers. If a node occurs on a path,
          then the time for which the node sequence persisted before another
          mutation is listed (for the first occurrence of the node if it
          occurs multiple times on the same path). Note that it *translateseqs*
          is True, then these are the persistence times to the first
          non-synonymous mutation, as only those mutations change the node
          sequences. The total length of the list for each node will be equal
          to the number of paths that contained that node.

        * *t.timetofirstmut* : if *translateseqs* is False, then this 
          is just equal to *t.persistence*. But if *translateseqs* is True,
          then the entries give the time after the occurrence of a node to 
          the first mutation of any time -- synonymous or nonsynonymous. In
          this case, entries in *t.timetofirstmut* will often be less than
          those in *t.persistence*, since the first mutation to a node will
          often by synonymous, which will change the nucleotide sequence
          but not the actual protein sequence node identity.

    To create a *Trajectory* object *t*, use the command::

        t = Trajectory(pathfile, translateseqs=False, printprogress=False)

    *pathfile* should be a string giving a name of an input file specifying one or
    more mutational paths. These files are of the format created by 
    ``mutpath_get_paths.py``. They must be readable by the *IteratePaths*
    function.

    *translateseqs* is an optional argument that is *False* by default. If it
    is set to *True*, then the sequences contained within *mutpathfile* are
    taken to represent coding nucleotide sequences, but the trajectory is
    built through protein sequence space. In other words, the nucleotide sequences
    in the paths are translated, and the trajectory is built from these
    translated sequences. All of the nodes and edges will therefore connect
    protein sequences. Note that no checking is made to ensure that the
    sequences translate properly: any stop codons are simply translated to '*',
    codons containing gaps are translated to '-', and sequences that do not have
    lengths that are multiples of three have the extra one or two nucleotides
    truncated.

    *printprogress* is a switch specifying that we print progress as we
    process paths. You may want to use this if you are processing a large
    number of paths and want to output the progress. By default it is *False*,
    meaning that nothing is printed. If you set it to an integer, it will then
    print to *sys.stdout* after processing every *printprogress* paths.
    """

    def __init__(self, pathfile, translateseqs=False, printprogress=False):
        """Intializes and returns a *Trajectory* object.

        Returns a *Trajectory* object *t* constructed from the collection of mutational
        paths encoded in *pathfile*. Detailed in main docstring for this class.
        """
        if not os.path.isfile(pathfile):
            raise IOError("Cannot find pathfile %s" % pathfile)
        self.npaths = 0
        self.nodes = {}
        self.edges = {}
        self.mutations = {}
        self.mutationstoca = {}
        self.persistence = {}
        if translateseqs:
            self.timetofirstmut = {}
        else:
            self.timetofirstmut = self.persistence
        self.startseq = self.endseq = None
        for (startseq, starttime, endseq, endtime, caseq, catime, tocamuts, fromcamuts) in IteratePaths(pathfile):
            onthispath = {}
            persistenceonthispath = {}
            timetofirstmutonthispath = {}
            self.npaths += 1
            if printprogress:
                if not (self.npaths % printprogress):
                    sys.stdout.write("Processed %d paths...\n" % self.npaths)
                    sys.stdout.flush()
            currentseq = list(startseq)
            nodetime = starttime
            if translateseqs:
                startseq = sequtils.Translate([('head', startseq)], readthrough_n=True, readthrough_stop=True, truncate_incomplete=True, translate_gaps=True)[0][1]
                endseq = sequtils.Translate([('head', endseq)], readthrough_n=True, readthrough_stop=True, truncate_incomplete=True, translate_gaps=True)[0][1]
            if self.startseq == None:
                self.startseq = startseq
                self.endseq = endseq
            assert self.startseq == startseq and self.endseq == endseq
            onthispath[startseq] = True
            if startseq in self.nodes:
                self.nodes[startseq] += 1
            else:
                self.nodes[startseq] = 1
            firstfromca = True
            for (mutlist, toca) in [(tocamuts, True), (fromcamuts, False)]:
                for (mut, time) in mutlist:
                    (wt, i, m) = (mut[0], int(mut[1 : -1]), mut[-1])
                    if not (1 <= i <= len(currentseq)):
                        raise ValueError("Position %d is out of range." % i)
                    if currentseq[i - 1] != wt:
                        raise ValueError("Identity mismatch for %s" % mut)
                    if wt == m:
                        raise ValueError("Invalid mutation %s" % mut)
                    s1 = ''.join(currentseq)
                    currentseq[i - 1] = m
                    s2 = ''.join(currentseq)
                    if translateseqs:
                        s1 = sequtils.Translate([('head', s1)], readthrough_n=True, readthrough_stop=True, truncate_incomplete=True, translate_gaps=True)[0][1]
                        s2 = sequtils.Translate([('head', s2)], readthrough_n=True, readthrough_stop=True, truncate_incomplete=True, translate_gaps=True)[0][1]
                    if not s2 in onthispath:
                        onthispath[s2] = True
                        if s2 in self.nodes:
                            self.nodes[s2] += 1
                        else:
                            self.nodes[s2] = 1
                    assert len(s1) == len(s2) == len(self.startseq) == len(self.endseq)
                    if self.persistence != self.timetofirstmut:
                        assert translateseqs
                        if s1 not in timetofirstmutonthispath:
                            timetofirstmutonthispath[s1] = True
                            if toca:
                                dt = time - nodetime
                            elif firstfromca:
                                dt = catime - nodetime + catime - time
                            else:
                                dt = nodetime - time
                            if s1 in self.timetofirstmut:
                                self.timetofirstmut[s1].append(dt)
                            else:
                                self.timetofirstmut[s1] = [dt]
                    if s1 != s2:
                        if s1 not in persistenceonthispath:
                            persistenceonthispath[s1] = True
                            if toca:
                                dt = time - nodetime
                            elif firstfromca:
                                firstfromca = False
                                dt = catime - nodetime + catime - time
                            else:
                                dt = nodetime - time
                            if s1 in self.persistence:
                                self.persistence[s1].append(dt)
                            else:
                                self.persistence[s1] = [dt]
                            nodetime = time
                        if translateseqs:
                            diffs = [i for i in range(len(s1)) if s1[i] != s2[i]]
                            assert len(diffs) == 1, str(diffs)
                            i = diffs[0]
                            mutstring = "%s%d%s" % (s1[i], i + 1, s2[i])
                        else:
                            mutstring = mut
                        if mutstring not in onthispath:
                            if mutstring in self.mutations:
                                self.mutations[mutstring].append(time)
                                if toca:
                                    self.mutationstoca[mutstring] += 1
                            else:
                                self.mutations[mutstring] = [time]
                                if toca:
                                    self.mutationstoca[mutstring] = 1
                                else:
                                    self.mutationstoca[mutstring] = 0
                            onthispath[mutstring] = True
                        tup = (s1, s2)
                        if not tup in onthispath:
                            onthispath[tup] = True
                            if tup in self.edges:
                                self.edges[tup] += 1
                            else:
                                self.edges[tup] = 1
            # check that path finished correctly
            if translateseqs:
                if sequtils.Translate([('head', ''.join(currentseq))], readthrough_n=True, readthrough_stop=True, truncate_incomplete=True, translate_gaps=True)[0][1] != endseq:
                    raise ValueError("Failed to end on endseq")
            elif ''.join(currentseq) != endseq:
                raise ValueError("Failed to end on endseq")
        if not self.npaths:
            raise ValueError("Failed to find any paths in %s" % pathfile)
        for key in self.nodes.iterkeys():
            if key != self.endseq:
                if len(self.persistence[key]) != self.nodes[key]:
                    raise ValueError("Incorect number of persistence entries")
            self.nodes[key] /= float(self.npaths)
        for key in self.edges.iterkeys():
            self.edges[key] /= float(self.npaths)
        if self.nodes[self.startseq] != 1:
            raise ValueError("The weight of startseq is not one")
        if self.nodes[self.endseq] != 1:
            raise ValueError("The weight of endseq is not one")



# Test with doctest
if __name__ == '__main__':
    import doctest
    doctest.testmod()

