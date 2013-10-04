"""Module for various functions for parsing phylogenetic trees.

This method is designed to take as input the text representations of
phylogenetic trees produced by programs such as BEAST, and convert them
to a form represented by the *Tree* and *Node* objects defined in the ``tree``
module.

Written by Jesse Bloom.


List of functions defined in this module
------------------------------------------

`TraceMutationPath` : Finds the mutational path between two sequences on a tree.

`AssignMutations` : Assigns mutations to node branches.

`BuildSeqsFromMutsAndTips` : Builds node sequences from tip and mutation information.

`RecursivelyAssignSeqsFromMuts` : assigns sequences using parental sequence and mutations.
`GetTreeString` : Gets the newick tree from a BEAST ``.tree`` file line.

`ReadTreeTranslateBlock` : Reads name translation block from ``.tree`` file.


Detailed documentation for functions
--------------------------------------
Detailed documentation for the functions is provided in their individual
docstrings below.

"""


import re
import tree
import sequtils


def TraceMutationPath(t, startstrain, endstrain):
    """Finds the mutational path between two sequences along a phylogenetic tree.

    Given a phylogenetic tree *t* with mutations assigned to branches and
    sequences assigned to nodes, this function traces the path along the
    tree from the tip node *startstrain* to the tip node *endstrain* and
    reconstructs the sequences along the way. This mutational path through
    sequence space is written to *pathfile*.

    CALLING VARIABLES:

    * *t* is a phylogenetic tree as a *tree.Tree* object. 
      All nodes should have mutations assigned to their branches in their
      *node.ancestorbranchinfo*, *node.rightdescendentinfo*, *node.leftdescendentinfo*
      attributes under the key 'mutations'. In addition, all nodes should have
      sequences assigned in their *node.info* under the key 'sequence'. In addition,
      all nodes should have their time preceding the most recent tip assigned
      in their *node.info* under the key 'time'. In practice,
      you can ensure this is true by performing the following two function
      calls before calling this function::
        
            parse_tree.AssignMutations(t.GetRoot())
            parse_tree.BuildSeqsFroMutsAndTips(t, allowambiguousmismatches)
            tree.AssignTimes(t)

      In addition, the tip nodes should have names assigned in their *node.info*
      dictionaries under the key 'strain'. These names are used to locate
      the tip nodes corresponding to *startstrain* and *endstrain*.

    * *startstrain* and *endstrain* are strings giving the names of the tip
      nodes that form the starting and ending sequences for the mutational path.
      There must be a unique tip node that has a *node.info* key 'strain' giving
      *startstrain* as its value -- this is the starting sequences. The same
      is true for *endstrain*.

    RESULT OF CALLING THIS FUNCTION

    This function returns a string specifying the mutational path and the timing of
    mutations. All times are measured in units before the most recent tip node
    (so the root node would have the largest value of time). Here is an example
    of a path that might be produced::

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

    In this output string, the first three entries give the name, sequence, and
    time for *startstrain*. The next three do the same for *endstrain*. The next
    two entries give the sequence and time of the most recent common ancestor
    of *startstrain* and *endstrain*. The mutational path is then specified,
    first tracing from *startstrain* to the common ancestor (so backwards along the
    tree), and then tracing from the common ancestor the *endstrain*. In each of
    these sections of the path, each line specifies a different mutation in 1, 2, ...
    numbering, followed by a colon and the time of the mutation (again measured
    backwards from the most recent tip). 
    """
    startnode = tree.GetNode(t, 'strain', startstrain)
    if not startnode:
        raise ValueError("Failed to find startstrain node %s" % startstrain)
    endnode = tree.GetNode(t, 'strain', endstrain)
    if not endnode:
        raise ValueError("Failed to find endstrain node %s" % startstrain)
    (commonancestor, pathlist) = tree.GetPath(startnode, endnode)
    assert startnode == pathlist[0]
    assert endnode == pathlist[-1]
    seq = list(startnode.info['sequence']) # keep track of sequence
    path_info = [
                 'startstrain_name %s' % startstrain,
                 'startstrain_seq %s' % startnode.info['sequence'],
                 'startstrain_time %f' % startnode.info['time'],
                 'endstrain_name %s' % endstrain,
                 'endstrain_seq %s' % endnode.info['sequence'],
                 'endstrain_time %f' % endnode.info['time'],
                 'commonancestor_seq %s' % commonancestor.info['sequence'],
                 'commonancestor_time %f' % commonancestor.info['time'],
                 'startstrain_to_commonancestor_path',
                ]
    movingforward = False # are we moving forward in the tree?
    reversemutations = []
    for node in pathlist:
        if node == commonancestor:
            reversemutations.sort() # sort to reverse chronological order
            reversemutations = ['%s%d%s : %f' % (m, i, wt, t) for (t, wt, i, m) in reversemutations]
            for mstring in reversemutations: # check sequence
                mstring = mstring.split()[0]
                (wt, i, m) = (mstring[0], int(mstring[1 : -1]), mstring[-1])
                if seq[i - 1] != wt:
                    raise ValueError("Tracing back to common ancestor, sequence mismatch for %s, actual identity is %s" % (mstring, seq[i - 1]))
                seq[i - 1] = m
            path_info += reversemutations
            path_info.append('commonancestor_to_endstrain_path')
            if movingforward:
                raise ValueError("We are moving forward, but just found the common ancestor. This shouldn't happen...")
            movingforward = True
        elif 'mutations' in node.ancestorbranchinfo:
            for (t, wt, i, m) in node.ancestorbranchinfo['mutations']:
                if movingforward:
                    path_info.append('%s%d%s : %f' % (wt, i, m, t))
                    if seq[i - 1] != wt:
                        raise ValueError("Tracing forward from common ancestor, sequence mismatch for %s%d%s, actual identity is %s" % (wt, i, m, seq[i - 1]))
                    seq[i - 1] = m
                else:
                    reversemutations.append((t, wt, i, m))
    if not movingforward:
        raise ValueError("We never moved forward from the common ancestor in the path trace. In general, this should not happen unless you assigned an endseq that is a direct ancestor of the starting sequence...")
    return '\n'.join(path_info)



def BuildSeqsFromMutsAndTips(t, allowambiguousmismatches):
    """Builds node sequences from tip and mutation information.

    This function takes a single argument *t* specifying a *tree.Tree*
    object. In order for this method to work, all nodes must have
    mutations assigned to their *branchinfo* dictionaries under the 
    key 'mutations' in the format given by **AssignMutations**. In other
    words, you should have already called::

        AssignMutations(t.GetRoot())

    In addition, all of the tip nodes of the tree should have their sequences
    (does not matter if these are protein or nucleic acid) specified in
    their *node.info* dictionaries under the key 'sequences'.

    This method uses the tip sequences and the mutation information to
    determine the sequences for all nodes of the tree, and places
    these in the *node.info* dictionary under the key 'sequence'. Furthermore,
    it makes sure that the tip sequences and mutations specified in the tree
    already are consistent, and raises an exception if they are not.

    First, it chooses a tip node and uses that to reconstruct the root node
    by tracing back up the tree to the root.
    Then it starts from this root node and recursively reconstruct all interior
    nodes and to check all tip nodes.

    *allowambiguousmismatches* is a list or set specifying any ambiguous characters
    that we might allow in tip nodes. In general, if a tip node contains an 
    ambiguous character, then during the sequence building a mismatch
    gaps as ambiguous characters, then during the sequence building a mismatch
    can occur. Any character that is present in *allowambiguousmismatches* does not
    count as a mismatch. For protein sequences, you might want:::
   
        allowambiguousmismatches = ['-', 'X']

    to allow ambiguous amino acids and gaps. For nucleotide sequences, you might
    want::

        allowambiguousmismatches = ['-', 'N']

    to allow ambiguous nucleotides and gaps.
    """
    # first get a tip node with as few ambiguous characters as possible
    (tip_list, internal_list) = ([], [])
    tree.ListsOfNodes(t.GetRoot(), tip_list, internal_list)
    if not tip_list:
        raise ValueError("no tip nodes")
    lowestsofar = None
    for tip in tip_list:
        try:
            seq = tip.info['sequence']
        except KeyError:
            raise KeyError("Failed to find 'sequence' key in node.info for tip node -- did you remember to build the tree so that tipnode sequences were assigned?")
        nambiguous = 0
        for ambiguous in allowambiguousmismatches:
            nambiguous += seq.count(ambiguous)
        if nambiguous == 0:
            (node, s) = (tip, list(seq))
            break
        elif lowestsofar == None or nambiguous < lowestsofar:
            (lowestsofar, node, s) = (nambiguous, tip, list(seq))
    else:
        raise ValueError("Failed to find a tip without ambiguous sequence")
    # now trace back from node to the root, reconstructing sequences
    while not node.root:
        muts = node.ancestorbranchinfo['mutations']
        # iterated in reversed order since muts is chronological
        for (t, wt, i, m) in reversed(muts):
            if s[i - 1] != m:
                if s[i - 1] not in allowambiguousmismatches:
                    raise ValueError("Mismatch trying to revert %s%d%s from %s" % (wt, i, m, s[i - 1]))
            s[i - 1] = wt
        node = node.ancestor
    assert node.root
    node.info['sequence'] = ''.join(s) # assign root sequence
    # now recursively build from the root to assign sequences
    RecursivelyAssignSeqsFromMuts(node, allowambiguousmismatches)



def RecursivelyAssignSeqsFromMuts(node, allowambiguousmismatches):
    """Assigns sequences using parental sequence and specified mutations.

    *node* is a *tree.Node* object.  *node* must have the following keys
    in its attribute dictionaries:

        * *node.info* must have the key 'sequence' and a value specifying
          a nucleotide or protein sequence for the node.

        * Unless *node* is a tip node (which is the case if *node.tip == True*),
          then *node.rightbranchinfo* and *node.leftbranchinfo* must both have the
          key 'mutations' specifying mutations in the format written by 
          **AssignMutations**. These mutations should be numbered in 1, 2, ...
          numbering of the sequence specified under 'sequence' in the *node.info*
          attribute.
   
    If *node* is a tip node (*node.tip == True*), this function does nothing.
    Otherwise, it uses the sequence specified in *node.info* and the mutations
    specified in *node.rightbranchinfo* and *node.leftbranchinfo* to build
    the sequences for the right and left descendents of *node*, and places
    these sequences under the key 'sequences' in the *info* attributes
    of these descendents. If a 'sequences' key already specifies a sequence
    for these descendents, then checks to make sure the sequences built from
    the ancestor and mutations matches this sequence -- if they do not,
    raises an exception. This function then recursively calls itself on the
    descendent nodes.

    If you call this function with *node* as the root node of a tree, this
    method will build the sequences for all nodes of the tree.

    *allowambiguousmismatches* has the same meaning as in **BuildSeqsFromMutsAndTips**.
    """
    if node.tip:
        return
    seq = node.info['sequence']
    for (descendent, muts) in [\
            (node.rightdescendent, node.rightbranchinfo['mutations']),\
            (node.leftdescendent, node.leftbranchinfo['mutations'])]:
        descendentseq = seq
        j = 0
        while j < len(muts):
            (t, wt, i, m) = muts[j]
            if descendentseq[i - 1] != wt:
                # I have sometimes seen the two mutations SIMULTANEOUS which
                # can lead to a mismatch here. This is how we deal with that...
                if (j + 1 < len(muts)) and t == muts[j + 1][0]:
                    # this mutation and the next are simultaneous, swap order
                    muts = muts[ : j] + [muts[j + 1]] + [muts[j]] + muts[j + 2 : ]
                    continue
                else:
                    raise ValueError("sequence mismatch, trying to make %s%d%s to %s (muts is %s)" % (wt, i, m, descendentseq[i - 1], str(muts)))
            if i == len(descendentseq):
                descendentseq = "%s%s" % (descendentseq[ : i - 1], m)
            else:
                descendentseq = "%s%s%s" % (descendentseq[ : i - 1], m, descendentseq[i : ])
            j += 1
        if 'sequence' in descendent.info:
            if descendent.info['sequence'] != descendentseq:
                if len(descendentseq) != len(descendent.info['sequence']):
                    raise ValueError('sequence length mismatch')
                mismatches = []
                for i in range(len(descendentseq)):
                    (x, y) = (descendentseq[i], descendent.info['sequence'][i])
                    if x != y:
                        if (y not in allowambiguousmismatches) and (x not in allowambiguousmismatches):
                            mismatches.append("%s%d%s" % (x, i + 1, y))
                if mismatches:
                    raise ValueError("Sequence mismatch:\n>descendentseq\n%s\n>descendent.info['sequence']\n%s\n\nMismatches:\n%s" % (descendentseq, descendent.info['sequence'], ', '.join(mismatches)))
        else:
            descendent.info['sequence'] = descendentseq
        RecursivelyAssignSeqsFromMuts(descendent, allowambiguousmismatches)




def AssignMutations(node):
    """Assigns mutations to node branches.

    CALLING VARIABLES

      * *nodes is a *tree.Node* object. Typically, you would call this function
        with *node* as the root node of a tree. The mutations should be specified
        in the *rightbranchinfo* and *leftbranchinfo* attributes of the nodes
        using the key 'history_all' as detailed below.

    RESULT OF THIS FUNCTION

    This function recursively proceeds down the tree descended from *node*
    and assigns mutations to the branches. 
    
    The *rightbranchinfo* and *leftbranchinfo* for each node are assumed
    to hold the information about the mutations. Specifically, if
    these dictionaries have an entry keyed by the string 'history_all',
    then this history should list all mutations (the absence of such
    a key indicates no mutations). These should be listed in forms
    such as the following::

        node.rightbranchinfo['history_all'] = '{{286,42.2151,A,S},{353,40.2638,S,C}}'

    where the first number is the position (residue or nucleotide number) in
    1, 2, ... numbering; the second number is the time of the mutation given
    in units such that the most current tip node has a value of 0 and the
    tree root has a value equal to the tree height, and the next numbers
    are the wildtype and then the mutant amino acids or nucleotides. This
    function creates an entry in the *branchinfo* dictionary that gives the
    following list, with the entries sorted from largest to smallest, which
    means that the list gives the mutations in chronological order since
    the largest times correspond to those in the most distant past::

        node.rightbranchinfo['mutations'] = [(42.2151, 'A', 286, 'S'), (40.2638, 'S', 353, 'C')]

    where specifically each entry in the list is the time since the root node, the 
    wildtype identity, the position, and the mutant identity. If there are no
    mutations, an empty list is added for 'mutations'.
    """
    if node.tip:
        return # nothing to do for tip nodes
    # handle descendents
    for (descendent, branchlength, branchinfo) in [(node.rightdescendent, node.rightbranch, node.rightbranchinfo), (node.leftdescendent, node.leftbranch, node.leftbranchinfo)]:
        mutations = []
        if 'history_all' in branchinfo:
            mutbrackets = branchinfo['history_all']
            assert mutbrackets[0] == '{' and mutbrackets[-1] == '}', mutbrackets
            mutbrackets = mutbrackets[1 : -1]
            if mutbrackets:
                assert mutbrackets[0] == '{' and mutbrackets[-1] == '}', mutbrackets
                for x in mutbrackets[1 : -1].split('},{'):
                    tup = x.split(',')
                    mutations.append((float(tup[1]), tup[2], int(tup[0]), tup[3]))
        mutations.sort()
        mutations.reverse()
        branchinfo['mutations'] = mutations
        AssignMutations(descendent)


def GetTreeString(line, replace_with_null=[]):
    """Gets the newick tree from a line from a BEAST ``.tree`` file.

    *line* is a string giving a line representing a tree from a BEAST
    tree file.

    *replace_with_null* is a list of regular expression patterns that we
    wish to remove from the tree string.  The replacements
    are done in the order that the items are specified. They are replaced
    with '' (nothing).

    Returns the newick tree portion of line only, with all patterns
    matching those listed in *replace_with_null* having been removed.
    """
    statematch = re.compile('^tree [A-Z]+\_{0,1}\d+ ')
    m = statematch.search(line)
    if not m:
        raise ValueError("Failed to find expected beginning of tree line:\n%s" % line)
    line = line[len(m.group(0)) : ]
    if line[ : 2] == '[&':
        i = tree.GetBalancedIndex(line, 0)
        line = line[i + 1 : ]
    if ' = [&R] ' == line[ : 8]:
        line = line[8 : ]
    elif '= [&R] ' == line[ : 7]:
        line = line[7 : ]
    else:
        raise ValueError("Failed to find expected information:\n%s" % line)
    for rematch in replace_with_null:
        line = re.sub(rematch, '', line)
    return line


def ReadTreeTranslateBlock(treefile):
    """Reads name translation block from ``.tree`` file.

    *treefile* is a string giving the name of a BEAST ``.trees`` file.

    A ``.trees`` file (such as that created by BEAST) contains a block
    called 'Translate' under the 'Begin trees' heading, such as::

        Begin trees;
            Translate
                1 'A/Aichi/2/1968_1968.50',
                2 'A/Akita/1/94_1994.50',
                3 'A/nt/60/1968_1968.50'
                ;

    This block is assigning the numbers used in the Newick strings to the
    full strain names. 

    This function returns a dictionary keyed by the integer strain codes
    and with the values being the strings giving the full name. The integer
    strain codes are given as strings, not numbers. For example
    the above block would return::

        {'1':'A/Aichi/2/1968_1968.50', '2':'A/Akita/1/94_1994.50', '3':'A/nt/60/1968_1968.50'}
    """
    d = {}
    begintrees = begintranslate = False
    linematch = re.compile("^(?P<num>\d+) \'(?P<name>.+)\'\,{0,1}$")
    for line in open(treefile):
        line = line.strip()
        if begintrees:
            if begintranslate:
                m = linematch.search(line)
                if not m:
                    if line == ';':
                        return d
                    else:
                        raise IOError("Error matching Translate block in %s" % treefile)
                else:
                    (num, name) = (m.group('num'), m.group('name'))
                    if num in d:
                        raise ValueError("Duplicate code of %d in %s Translate" % (num, treefile))
                    d[num] = name
            elif 'Translate' == line:
                begintranslate = True
        elif 'Begin trees;' == line:
            begintrees = True
    raise IOError("Failed to properly match Translate block in %s" % treefile)



# test with doctest
if __name__ == '__main__':
    import doctest
    doctest.testmod()
