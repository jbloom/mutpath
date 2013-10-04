"""Module establishing classes and functions for representing phylogenetic trees.

Written by Jesse Bloom.


Functions defined in this module
------------------------------------

`ApplyToNodes` : Applies a function recursively to all descendents of a node.

`GetBalancedIndex` : Finds the index of balancing parentheses character.

`GetPath` : Returns the path along a tree between two nodes.

`GetInfo` : Gets information from a BEAST format information bracket.

`TallySiteChangesAlongTree` : Counts times a site has changed along a tree.

`SplitNewick` : Splits a bifurcating tree into its next two branches.

`RecursivelyAssignDescendents` : Assigns descendents to a node of a tree.

`RecursivelySetNumbers` : Sets the numbers for nodes of a tree.

`RecursivelySetTipAndInternalNumbers` : Sets numbers for tip and internal nodes.

`GetNode` : Gets a specific node from a tree.

`ListsOfNodes` : Returns lists of all internal and tip nodes in a tree.

`QueryTree` : Gets information about a particular node of a tree.

`AssignTimes` : Assigns times before the most recent tip to all nodes of a tree.

`RecursivelyAssignTimes` : Assigns times to all descendents of a node.

`PrintAllTreeInformation` : Prints information about nodes, names, sequences.


Classes defined in this module
------------------------------------

`Tree` : A class for representing a phylogenetic tree (composed of *Node* objects).

`Node` : A class for representing the node objects that compose a *Tree*.


C-extensions for this module
----------------------------------
Some of the functions in this module are implemented in faster C code in the
``mutpath.ctree`` module. If this module is available, by default these faster implementations
are used. The C implementations are automatically called when calling the Python 
function, so there is no need to call functions in ``mutpath.ctree`` directly.

The following functions can utilize fast C extensions:

    * **GetBalancedIndex**


Details of functions and classes
----------------------------------
Detailed descriptions of the functions and classes are given
in their individual docstrings below.


"""



import sys
import re
import warnings
try:
    import mutpath.ctree
    _ctreeimported = True
except ImportError:
    warnings.warn("Cannot import mutpath.ctree. Will use slower pure Python implementations.")
    _ctreeimported = False


def ApplyToNodes(node, F):
    """Applies a function recursively to all descendents of a node.

    *node* is a `Node` object.

    *F* is a function that takes as input a node.  This function applies
    *F* to *node* and to all of its descendents.

    Typically, you will call this function with *node* as the root
    node of a tree, and then the function will be applied to the 
    whole tree.
    """
    F(node)
    if not node.tip:
        ApplyToNodes(node.rightdescendent, F)
        ApplyToNodes(node.leftdescendent, F)


def GetBalancedIndex(s, i, use_ctree=True):
    """Finds the index of the balancing parentheses character.

    *s* is a string.

    *i* is an integer index (*0 <= i < len(s)*).  *s[i]* must be one of the characters
    "(", ")", "[", "]", "{", or "}" -- otherwise an exception is raised.

    *use_ctree* specifies that we use the faster implementations in the ``ctree``
    C extension if available. Is *True* by default; if you set to *False* will
    use the pure Python code.

    If *s[i]* is start bracket character such as "(", then moves
    forward in the string to find and return the index of the
    balancing bracket character.  If *s[i]* is an end bracket character,
    then moves backward in the string to find and return the index
    of the balancing bracket character.  Raises an exception of the
    character cannot be balanced.
    
    Examples:

    >>> GetBalancedIndex('hi (what did you say)!', 3)
    20

    >>> GetBalancedIndex('hi (what did you say)!', 20)
    3

    >>> GetBalancedIndex('(507:12.32[&h="A,(B)"]),(234)', 0)
    22

    >>> GetBalancedIndex('(507:12.32[&h="A,(B)"]),(234)', 10)
    21
    
    """
    if use_ctree and _ctreeimported:
        return mutpath.ctree.GetBalancedIndex(s, i)
    if not (0 <= i < len(s)):
        raise ValueError("Invalid index of %d for string of length %d" % (i, len(s)))
    if s[i] in ['(', '[', '{']:
        partner = {'(':')', '[':']', '{':'}'}[s[i]]
        j = i + 1
        nested = 0
        while j < len(s):
            if s[j] == s[i]:
                nested += 1
            elif s[j] == partner:
                if nested > 0:
                    nested -= 1
                else:
                    return j
            j += 1
        raise ValueError("Failed to balance %s" % s[i])
    elif s[i] in [')', ']', '}']:
        partner = {')':'(', ']':'[', '}':'{'}[s[i]]
        j = i -1
        nested = 0
        while j >= 0:
            if s[j] == s[i]:
                nested += 1
            elif s[j] == partner:
                if nested > 0:
                    nested -= 1
                else:
                    return j
            j -= 1
        raise ValueError("Failed to balance %s" % s[i])
    else:
        raise ValueError("Character is not a bracket: %s" % s[i])


def GetPath(startnode, endnode, markpath=None):
    """Returns the path along the tree between two nodes.

    *startnode* and *endnode* are two nodes from the same *Tree*
    object. This function traces the path from *startnode* to *endnode*. 

    The returned variable is 2-tuple of the form
    *(commonancestor, pathlist)*
    where *commonancestor' is the last common ancestor node, and
    *pathlist* is a list of the form
    *[startnode, node1, node2, ..., endnode]*
    where *node1*, *node2*, etc represent intermediate nodes
    between *startnode* and *node1*.  Initially, each
    node in the list will be the immediate ancestor of its predecessor
    in the list, until we have traced back to the last common ancstor.
    After that, each node in the list will be the descendent of its
    predecessor in the list.  

    *markpath* gives us the option to mark the nodes and branches
    by writing in their information dictionaries to indicate
    they are on the path.  If *markpath* is set to a non-null
    string, then for each node and branch on the path, an
    entry is created in the dictionary keyed by that
    string, and with the integer value of one.
    """
    if startnode == endnode:
        raise ValueError("startnode and endnode are equal.")
    startnode_ancestors = [startnode]
    while not startnode_ancestors[-1].root:
        startnode_ancestors.append(startnode_ancestors[-1].ancestor)
    endnode_ancestors = [endnode]
    while not endnode_ancestors[-1].root:
        endnode_ancestors.append(endnode_ancestors[-1].ancestor)
    endnode_ancestors.reverse()
    for i in range(len(startnode_ancestors)):
        if startnode_ancestors[i] in endnode_ancestors:
            break
    else:
        raise ValueError('Failed to find common ancestor.')
    commonancestor = startnode_ancestors[i]
    if markpath:
        assert isinstance(markpath, str)
        for node in startnode_ancestors[ : i]:
            node.info[markpath] = 1
            node.ancestorbranchinfo[markpath] = 1
        for node in endnode_ancestors[endnode_ancestors.index(commonancestor) + 1 : ]:
            node.info[markpath] = 1
            node.ancestorbranchinfo[markpath] = 1
        commonancestor.info[markpath] = 1
    pathlist = startnode_ancestors[ : i] + endnode_ancestors[endnode_ancestors.index(commonancestor) : ]
    return (commonancestor, pathlist)


def GetInfo(s):
    """Gets information from a information bracket.

    *s* is a string that defines a valid information bracket in the
    BEAST annotation format.  That is,
    it must start with '[&' and end with ']'.  These two brackets
    must represent balanced parentheses.  There can then be an
    arbitrary number of pieces of information defined, in the form
    *key=value* where *value* is either a number or a string.  

    Returns a dictionary keyed by the keys as strings, and with the values as
    strings.  If the value represents an integer, is converted into one. If the
    value represents a float, is converted into one. Otherwise, is kept as a
    string. If the string is enclosed in brackets, returns everything
    in the brackets.

    Currently will not handle keys or values that contain the "=" character,
    and will not handle values with unbalanced parentheses.
    Raises a *ValueError* if the string cannot be parsed in this fashion.
    Also will raise a *ValueError* if it finds duplicate keys.
    
    Example:
  
    >>> s = '[&history_1={A:G:2.33},seq="ATGC",d=10.23,i=9]'
    >>> GetInfo(s) == {'history_1':'{A:G:2.33}', 'seq':'ATGC', 'd':10.23, 'i':9}
    True
    """
    if not (isinstance(s, str) and len(s) >= 3 and s[ : 2] == '[&' and s[-1] == ']'):
        raise ValueError("Was not passed a valid information string: %s" % str(s))
    s = s[2 : -1]
    info = {}
    keymatch = re.compile("^\S+?\=")
    m = keymatch.search(s)
    intmatch = re.compile("^\-{0,1}\d+($|\,)")
    floatmatch = re.compile("^\-{0,1}\d+\.\d*(E\-\d+){0,1}($|\,)")
    stringmatch1 = re.compile("^'[^']*'($|\,)")
    stringmatch2 = re.compile('^"[^"]*"($|\,)')
    while m:
        key = m.group(0)
        s = s[s.index(key) + len(key) : ]
        if key[0] == ',':
            if len(key) <= 2:
                raise ValueError("Problem with key:\n%s" % key)
            key = key[1 : -1]
        else:
            key = key[ : -1]
        if intmatch.search(s):
            svalue = intmatch.search(s).group(0)
            if svalue[-1] == ',':
                value = int(svalue[ : -1])
            else:
                value = int(svalue)
        elif floatmatch.search(s):
            svalue = floatmatch.search(s).group(0)
            if svalue[-1] == ',':
                value = float(svalue[ : -1])
            else:
                value = float(svalue)
        elif stringmatch1.search(s):
            svalue = stringmatch1.search(s).group(0)
            if svalue[-1] == ',':
                value = svalue[1 : -2]
            else:
                value = svalue[1 : -1]
        elif stringmatch2.search(s):
            svalue = stringmatch2.search(s).group(0)
            if svalue[-1] == ',':
                value = svalue[1 : -2]
            else:
                value = svalue[1 : -1]
        elif s[0] in ['[', '{', '(']:
            j = GetBalancedIndex(s, 0)
            if j + 1 < len(s) and s[j + 1] == ',':
                svalue = value = s[ : j + 1]
            elif j + 1 == len(s):
                svalue = value = s[ : j + 1]
            else:
                raise ValueError("Problem with:\n%s" % s)
        else:
            raise ValueError("Cannot parse value from:\n%s" % s)
        if key in info:
            raise ValueError("Duplicate keys of %s" % key)
        else:
            info[key] = value
        s = s[len(svalue) : ]
        m = keymatch.search(s)
    if s:
        raise ValueError("Unparsed string portion: %s" % s)
    return info


def TallySiteChangesAlongTree(t):
    """Counts the number of times a site has changed along a reconstructed tree.

    Takes as input a single argument *t* which shuld represent a *Tree*
    object.  This tree must be rooted, and must have sequence assigned to
    both all tip nodes and all interior nodes.  That is, the interior
    sequences must already be reconstructed by some method.  These sequences
    should all be the same length (i.e., they should be aligned).

    This method traverses along the tree starting at the root.  For each site,
    it counts the number of times that the identity of that site has
    changed, and what these changes have been.  It returns a dictionary *tally*
    keyed by sequence numbers (ranging from 1 to the length of the sequence).
    Entry *tally[i]* is the 2-tuple *(n, changes)*.  *n* is the total number
    of times site *i* is inferred to have changed identities.  *changes* is
    itself a dictinary.  For each transition from *x* to *y* at site *i*, 
    *changes[(x, y)]* counts the number of times that transition occurred.
    """
    raise ValueError("Not updated since changing node sequences from node.seq to node.info['sequence']")
    assert isinstance(t, Tree)
    root = t.GetRoot()
    assert root.root
    if root.tip:
        raise ValueError("Root is also tip, indicating a null tree.")
    if not root.seq:
        raise ValueError("No sequence defined for root node.")
    seqlength = len(root.seq)
    tally = {}
    for i in range(1, seqlength + 1):
        tally[i] = (0, {})
    
    def _RecursivelyTally(tally, ancestorseq, node, seqlength):
        """Recursive internal subroutine for tallying changes."""
        if not node.seq:
            raise ValueError("Node has no sequence.")
        if len(node.seq) != seqlength:
            raise ValueError("Node sequence has different length.")
        for i in range(seqlength):
            if ancestorseq[i] != node.seq[i]:
                (n, itally) = tally[i + 1]
                n += 1
                tup = (ancestorseq[i], node.seq[i])
                if tup in itally:
                    itally[tup] += 1
                else:
                    itally[tup] = 1
                tally[i + 1] = (n, itally)
        if node.tip:
            return
        else:
            _RecursivelyTally(tally, node.seq, node.rightdescendent, seqlength)
            _RecursivelyTally(tally, node.seq, node.leftdescendent, seqlength)

    # Start recursively tallying down the tree
    _RecursivelyTally(tally, root.seq, root.rightdescendent, seqlength)
    _RecursivelyTally(tally, root.seq, root.leftdescendent, seqlength)
    return tally



def SplitNewick(newick_tree):
    """Splits a bifurcating phylogenetic tree into its next two branches.

    *newick_tree* is a string specifying a sub-tree of a Newick tree.
    If the string specifies the full tree, the trailing semicolon
    and the parentheses surrounding the root must already be
    removed.  In other words, the tree::

        (((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0.0;

    is not acceptable, but the tree::

        ((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7

    is acceptable.  The same limitations to the Newick tree apply
    as are listed in the documentation string for the *__init__*
    method of *Tree*.

    Returns the 2-tuple *(left_branch, right_branch)*.  The elements of this
    tuple specify the left and right branches of the tree, respectively.
    Each branch is itself a 5-tuple of the form 
    *(branch, length, tip, node_info, branch_info)*, where:

        * *branch* is a string giving the subtree if this branch is a subtree,
          or the name of the node if it is a tip.  

        * *length* is a number giving
          the length of the branch leading to this tip or subtree.
          
        * *tip* is *True* if the branch is a tip or *False* otherwise.

        * *node_info* gives any information specified for the node.
        
        * *branch_info* gives any information specified for the branch.

    Example:

    >>> newick_tree = '((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7'
    >>> (left_branch, right_branch) = SplitNewick(newick_tree)
    >>> right_branch == ("Five", 0.7, True, {}, {})
    True
    >>> left_branch == ("(One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2", 0.3, False, {}, {})
    True
    
    """
    floatmatch = re.compile('^\-{0,1}\d+\.{0,1}\d*(E\-\d+){0,1}')
# JDB, March-16-2013, replacing endfloatmatch with rfind for speed
#    endfloatmatch = re.compile('-{0,1}\d+\.{0,1}\d*(E\-\d+){0,1}$')
    assert isinstance(newick_tree, str)
    # We find the comma that is not enclosed in any parentheses.  There should
    # be exactly one such comma.
    if newick_tree[0] == '(':
        j = GetBalancedIndex(newick_tree, 0)
        icomma = j + 1
    elif newick_tree[-1] == ')':
        j = GetBalancedIndex(newick_tree, len(newick_tree) - 1)
        icomma = j - 1
    else:
        # move back from last
#        m = endfloatmatch.search(newick_tree)
#        if not m:
#            raise ValueError("Problem parsing the tree string - no parentheses, but does not end in branch length:\n%s" % newick_tree)
#        j = len(newick_tree) - 1 - len(m.group(0))
        j = max(newick_tree.rfind(']'), newick_tree.rfind(':'))
        if j == -1:
            raise ValueError("Problem parsing the tree string - no parentheses, but does not end in branch length:\n%s" % newick_tree)
        if newick_tree[j] == ']':
            # branch annotation
            j = GetBalancedIndex(newick_tree, j)
            assert newick_tree[j : j + 2] == '[&'
            j -= 1
        if newick_tree[j] != ':':
            raise ValueError("Problem parsing tree, expected colon:\n%s" % newick_tree)
        if newick_tree[j - 1] == ']':
            # node annotation
            j = GetBalancedIndex(newick_tree, j - 1)
            assert newick_tree[j : j + 2] == '[&'
        j -= 1 # JDB, moved out of preceding if clause on 3-16-13 to fix apparent bug
        if newick_tree[j] == ')':
            j = GetBalancedIndex(newick_tree, j)
        while newick_tree[j] != ',':
            j -= 1
            if j < 0:
                raise ValueError("Never found expected comma in:\n%s" % newick_tree)
        icomma = j
    if newick_tree[icomma] != ',':
        if newick_tree[icomma : icomma + 2] == '[&': # node information
            j = GetBalancedIndex(newick_tree, icomma)
            icomma = j + 1
        if newick_tree[icomma] != ':': # branch length
            raise ValueError("Expected to find branch length:\n%s" % newick_tree)
        icomma += 1
        if newick_tree[icomma : icomma + 2] == '[&': # branch information
            j = GetBalancedIndex(newick_tree, icomma)
            icomma = j + 1
        m = floatmatch.search(newick_tree[icomma : ])
        if m:
            icomma += len(m.group(0))
        else:
            raise ValueError("Failed to match branch length, listing icomma, newick_tree[icomma : ], and newick_tree:\n%d\n%s\n%s" % (icomma, newick_tree[icomma : ], newick_tree))
        if newick_tree[icomma] != ',':
            raise ValueError("Failed to find the expected comma:\n%s" % newick_tree)
    (branch_l, branch_r) = (newick_tree[ 0 : icomma], newick_tree[icomma + 1 : ])
    left_branch = []
    right_branch = []
    for (in_branch, out_branch) in [(branch_l, left_branch), (branch_r, right_branch)]:
#        m = endfloatmatch.search(in_branch)
#        if not m:
#            raise ValueError("Failed to find branch length at the end:\n%s" % in_branch)
#        j = len(in_branch) - 1 - len(m.group(0))
#        length = float(m.group(0))
        j = max(in_branch.rfind(']'), in_branch.rfind(':'))
        if j == -1:
            raise ValueError("Failed to find branch length at the end:\n%s" % in_branch)
        length = float(in_branch[j + 1 : ])
        if in_branch[j] == ']':
            # there is branch info here
            i = GetBalancedIndex(in_branch, j)
            assert in_branch[i : i + 2] == '[&'
            branch_info = GetInfo(in_branch[i : j + 1])
            j = i - 1
        else:
            branch_info = {}
        if in_branch[j] != ':':
            raise ValueError("Expected a colon in:\n%s\nInstead found:\n%s\nSeparating:\n%s\nAnd\n%s" % (in_branch, in_branch[j], in_branch[ : j], in_branch[j : ]))
        if in_branch[j - 1] == ']':
            # there is node info here
            i = GetBalancedIndex(in_branch, j - 1)
            assert in_branch[i : i + 2] == '[&'
            node_info = GetInfo(in_branch[i : j])
            j = i 
        else:
            node_info = {}
        branch = in_branch[ : j]
        if branch[0] == '(':
            assert branch[-1] == ')'
            branch = branch[1 : -1]
            tip = False
        else:
            tip = True
        out_branch += [branch, length, tip, node_info, branch_info]
    return (tuple(left_branch), tuple(right_branch))



def RecursivelyAssignDescendents(node, newick_tree, scalebranchlength, tipnames_dict, make_lengths_non_negative=False):
    """Recursively assigns the descendents to a node in a phylogenetic tree.

    This method assigns the descendent nodes for *node*, and does the same
    for these descendent nodes, until the tips of the tree have
    been reached.  The tips of the trees have names assigned.  Branch
    lengths are assigned to all nodes.

    CALLING VARIABLES:

    *node* is a *Node* that is not a tip and has no descendents assigned.
        
    *newick_tree* is a string giving the Newick tree
    that has *node* as its root.  The same restrictions apply
    to *newick_tree* as those mentioned in the *__init__* method
    of *Tree*.  

    *scalebranchlength* is an optional argument that allows us to rescale
    all branch lengths.  If it has a value of *None*, then no
    rescaling is done.  If it has another value, then it should be a number
    greater than zero.  In this case, all branch lengths specified by 
    *newick_tree* are multiplied by *scalebranchlength*.

    *tipnames_dict* has a similar meaning as described in the *__init__*
    method of *Tree*.  It is used to assign properties to the tip nodes.  It
    can still be *None*, meaning that no sequences are assigned.  Otherwise
    it is a dictionary, with the keys being the nodes names and the values
    themselves being dictionaries giving key / value pairs to add to the
    *node.info* dictionaries. There must be a key for every tip node.

    *makes_lengths_non_negative* is a Boolean switch specifying that we set
    any branch lengths that are less than zero to instead be zero.
    This is done if it is *True*.

    RESULT OF CALLING THIS FUNCTION:

    At the conclusion of this function, *node* will have its descendents
    assigned as specified by *newick_tree*. *node* will therefore represent
    the root of a tree.
    """
    if scalebranchlength == None:
        scalebranchlength = 1.0 # No scaling is the same as scaling by one
    else:
        assert 0 < scalebranchlength and isinstance(scalebranchlength, (int, float))
        pass # we scale by whatever number is set
    assert isinstance(node, Node) and not node.tip
    (left_branch, right_branch) = SplitNewick(newick_tree)
    if make_lengths_non_negative and left_branch[1] < 0:
        left_branch = (left_branch[0], 0.0, left_branch[2], left_branch[3], left_branch[4])
    if make_lengths_non_negative and right_branch[1] < 0:
        right_branch = (right_branch[0], 0.0, right_branch[2], right_branch[3], right_branch[4])
    # set left branch 
    if left_branch[2]: # left branch is a tip
        name = left_branch[0]
        info = left_branch[3]
        if tipnames_dict:
            if name not in tipnames_dict:
                raise KeyError("Failed to find %s in tipnames_dict" % name)
            for (key, value) in tipnames_dict[name].iteritems():
                if key in info:
                    raise KeyError("Duplicate key of %s" % key)
                info[key] = value
        left_node = Node(tip=True, ancestor=node, name=name, ancestorbranch=left_branch[1] * scalebranchlength, info=info, ancestorbranchinfo=left_branch[4])
    else: # left branch is not a tip, proceed recursively
        left_node = Node(ancestor=node, ancestorbranch=left_branch[1] * scalebranchlength, info=left_branch[3], ancestorbranchinfo=left_branch[4])
        RecursivelyAssignDescendents(left_node, left_branch[0], scalebranchlength, tipnames_dict, make_lengths_non_negative=make_lengths_non_negative)
    node.leftdescendent = left_node
    node.leftbranch = left_branch[1] * scalebranchlength
    node.leftbranchinfo = left_branch[4]
    # set right branch 
    if right_branch[2]: # left branch is a tip
        name = right_branch[0]
        info = right_branch[3]
        if tipnames_dict:
            if name not in tipnames_dict:
                raise KeyError("Failed to find %s in tipnames_dict" % name)
            for (key, value) in tipnames_dict[name].iteritems():
                if key in info:
                    raise KeyError("Duplicate key of %s" % key)
                info[key] = value
        right_node = Node(tip=True, ancestor=node, name=name, ancestorbranch=right_branch[1] * scalebranchlength, info=info, ancestorbranchinfo=right_branch[4])
    else: # right branch is not a tip, proceed recursively
        right_node = Node(ancestor=node, ancestorbranch=right_branch[1] * scalebranchlength, info=right_branch[3], ancestorbranchinfo=right_branch[4])
        RecursivelyAssignDescendents(right_node, right_branch[0], scalebranchlength, tipnames_dict, make_lengths_non_negative=make_lengths_non_negative)
    node.rightdescendent = right_node
    node.rightbranch = right_branch[1] * scalebranchlength
    node.rightbranchinfo=right_branch[4]



def RecursivelySetNumbers(node, number):
    """Function for setting the numbers for nodes of a tree.
   
    Calling this function on the root node of a tree will set unique 
    numbers for all nodes of the tree.

    A node's number is the value of *node.number*.

    Calling this function with *node* as the root node and *number* set to zero
    will return an integer giving the number of nodes in a tree (tip and
    internal nodes combined).  Typically, that is how you will use
    the function:: 
    
        RecursivelySetNumbers(tree.GetRoot(), 0)

    On call, *node* should be a *Node* object.

    *number* should be an integer that is the number that is set for *node*.

    Sets a number for *node*, and also recursively sets numbers for any of
    its descendents.  The returned value is the next number (number + 1)
    that is not set for any node.
    """
    node.number = number
    number += 1 # number of next node
    if not node.tip:
        number = RecursivelySetNumbers(node.rightdescendent, number)
        number = RecursivelySetNumbers(node.leftdescendent, number)
    return number



def RecursivelySetTipAndInternalNumbers(node, tip_number, internal_number):
    """Function for setting numbers for tip and internal nodes.

    The node numbers are the values accessed by *node.number*.  Any
    existing values for these numbers are cleared.

    Typically this function will initially be called with *node* as the
    root node of a tree, and *tip_number* and *internal_number* both
    equal to zero.  In this case, all of the tip nodes in the tree
    will be assigned integer numbers going from 0, 1, ...
    All of the inernal nodes in the tree will also be assigned integer
    numbers going from 0, 1, ...  Crucially, the numbers for the internal
    nodes are assigned so that if *nodej* is a descendent of *nodei*, then
    *nodej.number < nodei.number*. This leads to the implication
    that the root node will have the smallest internal node number.  Any existing 
    numbers for the nodes are cleared.  Note that it is NOT the case that all nodes
    in the tree have unique numbers; all tip nodes have unique numbers,
    and all internal nodes have unique numbers.  
    
    The returned value is
    the 2-tuple *(n_tips, n_internal)* where *n_tips* is the number of
    tip nodes and *n_internal* is the number of internal nodes.

    In its implementation, this function works by recursively calling itself.  In
    these recursive calls, the values of *tip_number* and *internal_number* are
    incremented as nodes are assigned numbers -- the value of each of these
    corresponds to the next number that will be assigned to the next encountered
    tip or internal node.
    """
    if node.tip:
        node.number = tip_number
        return (tip_number + 1, internal_number)
    else:
        (tip_number, internal_number) = RecursivelySetTipAndInternalNumbers(node.rightdescendent, tip_number, internal_number)
        (tip_number, internal_number) = RecursivelySetTipAndInternalNumbers(node.leftdescendent, tip_number, internal_number)
        node.number = internal_number
        return (tip_number, internal_number + 1)

def GetNode(t, key, value):
    """Gets a specific node from a tree.

    *t* is a *tree.Tree* object.

    *key* and *value* are used to specify how we locate the node. For nodes
    in the tree *t*, we look through the *node.info* attribute for each node
    for a key that matches *key* and specifies a value equal to *value*.
    This function returns the first node it finds where *node.info[key]*
    exists and specifies *value*. If no such node is found, returns None.
    """
    tip_list = []
    internal_list = []
    ListsOfNodes(t.GetRoot(), tip_list, internal_list)
    for node in tip_list + internal_list:
        if (key in node.info) and (node.info[key] == value):
            return node
    return None


def ListsOfNodes(node, tip_list, internal_list):
    """Function that returns lists of all tip and internal nodes in a tree.

    Typically you will call this function with *node* as a *Node* object giving the 
    root node of a tree, and *tip_list* and *internal_list* as set to empty lists.

    On return, *tip_list* is a list of all tip nodes in the tree with *node* as its 
    root node.  *internal_list* is a list of all internal nodes in the tree with
    *node* as its root.

    If the method is called with *node* equal to something other than the root node,
    then the lists are of all nodes in the subtree with root *node*.

    If the method is called with *tip_list* and *internal_list* as non-empty lists,
    then the tip and internal nodes are merely appended to the existing entries
    of these lists.

    This function works by recursively calling itself.
    """
    if node.tip:
        tip_list.append(node)
    else:
        internal_list.append(node)
        ListsOfNodes(node.rightdescendent, tip_list, internal_list)
        ListsOfNodes(node.leftdescendent, tip_list, internal_list)



def QueryTree(tree, nodelocation):
    """Gets information about a particular node of a tree.

    *tree* is a *Tree* object specifying a phylogenetic tree.

    *nodelocation* is a string giving the location of the node that
    we want to query from the tree.  It is specified in terms of
    tracing the tree forward from its root node.  A node location
    of '' (empty string) refers to the root node.  A node location
    of 'right' refers to the right descendent of the root node.
    A node location of 'right, left, left' refers to the node reached
    by beginning at the root node, taking the right descendent node,
    then this node's left descendent, then that node's left descendent.

    If *nodelocation* fails to specify an existing node of the tree, then
    an exception is raised.

    If *nodelocation* specifies an existing node, the function returns the 3-tuple:
    *(tip, name, sequence)* where:

        * *tip* is *True* if the node is a tip node, and *False* otherwise.

        * *name* is the string giving the name of the node, or *None* if
          no name is assigned to the node.

        * *sequence* is a string giving the sequence assigned to the node,
          or *None* if no sequence is assigned to this node.
    """
    raise ValueError("Not updated since node.seq moved to node.info['sequence']")
    assert isinstance(tree, Tree)
    assert isinstance(nodelocation, string)
    nodelocation = [entry.strip() for entry in nodelocation.split(',')]
    node = tree.GetRoot()
    for direction in nodelocation:
        if node.tip:
            raise ValueError("'nodelocation' is specifying a descendent to a tip node with a %s" % direction)
        if direction == 'right':
            node = node.rightdescendent
        elif direction == 'left':
            node = node.leftdescendent
        else:
            raise ValueError("Invalid descendent direction of %s.  Must be 'right' or 'left'." % direction)
    return (node.tip, node.name, node.sequence)


def AssignTimes(t):
    """Assigns times to all nodes of a tree.

    *t* is a *tree.Tree* object for which we are assigning node times.

    At the conclusion of this function, all nodes of *t* have times
    assigned in their *node.info* attributes under the key 'time'.
    The value for this key is a number that gives the time for that 
    node. Times are measured backwards from the most recent tip node.
    So the most recent tip node has a time of 0.0, and the root node
    has the largest value of time, which is the total height of the tree.

    In addition, all nodes will also have a key in their *node.info* attribute
    of 'time_since_root' which gives the time of the node since the root.
    The value of this key is 0.0 for the root node, and is the total height
    of the tree for the most recent tip node.
    """
    treeheight = RecursivelyAssignTimes(t.GetRoot(), 0.0, 'time_since_root')
    (tip_list, internal_list) = ([], [])
    ListsOfNodes(t.GetRoot(), tip_list, internal_list)
    for node in tip_list + internal_list:
        node.info['time'] = treeheight - node.info['time_since_root']



def RecursivelyAssignTimes(node, origin_time, time_label):
    """For all descendents of a node, assigns the time since this node.

    *node* is a *Node object*.  Assigns this node an entry in *node.info*
    of with the key *time_label* (a string) and the value *origin_time*   
    (a number). Then recursively assigns each descendent node a time equal
    to *origin_time* plus the length of the branch to that node.

    Typically, you might call this method with *node* equal to the
    root node of a tree, *origin_time* equal to 0.0, and *time_label*
    equal to "time_since_root".  Then::
    
        node.info["time_since_root"] = 0.0

    and all descendent nodes will have the time since the root specified.
    In addition, the function returns a single number representing
    the largest value of the time for any node in the subtree. This 
    number therefore represents the total height of the subtree rooted
    at node (if *origin_time* is zero on the call).
    """
    node.info[time_label] = origin_time
    if node.tip:
        return origin_time
    right_time = node.rightbranch + origin_time
    x = RecursivelyAssignTimes(node.rightdescendent, right_time, time_label)
    left_time = node.leftbranch + origin_time
    y = RecursivelyAssignTimes(node.leftdescendent, left_time, time_label)
    return max(x, y, origin_time)


def PrintAllTreeInformation(tree, f=sys.stdout):
    """Prints information about all nodes, names, and sequences in a tree.

    *tree* is a *Tree* object specifying a phylogenetic tree.

    *f* is a writeable filelike object.  By default, it is set to *sys.stdout*.
    You can also set it to some other writeable file.

    This function writes to *f* a FASTA file containing an entry for every node
    in *tree*.  The headers begin with the node's location in the tree,
    then whether this is an interior or tip node, and finally the node's
    assigned name.  The sequence is the value in *node.info['sequence']* for each
    node. Each node is assumed to have an entry with this key of 'sequence'
    in its *node.info* dictionary. The nodes location is written
    as "ROOT" if the node is a root node, and otherwise as a string of
    comma delimited "LEFT" and "RIGHT" characters indicating how the node
    is reached from the root node.  If no name is assigned for the node,
    the name is printed as "None". For example, the function might write::

        >ROOT INTERIOR None 
        None
        >LEFT INTERIOR None
        None
        >LEFT,LEFT TIP A/WSN/33 (H1N1)
        MTAGCKL
        >LEFT,RIGHT TIP A/PR/8/34 (H1N1)
        MTAGNL
        >RIGHT TIP A/Aichi/1968 (H3N2)
        MSGGNL
    """
    root = tree.GetRoot()
    assert not root.tip
    f.write(">ROOT INTERIOR %s\n%s\n" % (root.name, root.info['sequence']))
    # create function for recursively writing out node information
    def _RecursivelyWriteNode(node, location):
        """Recursive function for printing node descriptions."""
        if node.tip:
            f.write(">%s TIP %s\n%s\n" % (location, node.name, node.info['sequence']))
        else:
            f.write(">%s INTERIOR %s\n%s\n" % (location, node.name, node.info['sequence']))
            _RecursivelyWriteNode(node.leftdescendent, "%s,LEFT" % location)
            _RecursivelyWriteNode(node.rightdescendent, "%s,RIGHT" % location)
    # write out all of the descendent nodes
    _RecursivelyWriteNode(root.leftdescendent, 'LEFT')
    _RecursivelyWriteNode(root.rightdescendent, 'RIGHT')



class Tree(object):
    """Class for representing a phylogenetic tree.

    Currently only works for representing a rooted bifurcating tree.

    The tree is composed of *Node* objects.
    
    Calling arguments are documented in detail in the *__init__* method.
    """

    def __init__(self, newick_tree, scalebranchlength=None, tipnames_dict=None, make_lengths_non_negative=False):
        """Creates a phylogenetic tree.

        This class can currently represent only rotted bifurcating trees.

        CALLING VARIABLES: 
        
        *newick_tree* should be a string specifying a tree in Newick
        format with branch lengths.  The tree tips must have non-empty names,
        and the names cannot contain spaces and must be alphanumeric.  
        Names cannot be duplicated. Any spaces and line breaks are parsed out of
        *newick_tree*. Currently, only certain aspects of the full Newick format
        (http://evolution.gs.washington.edu/phylip/newick_doc.html) are
        handled by this method.  Specifically, there is no mechanism for:

            * handling names for internal nodes

            * handling quoted labels or labels with non alphanumeric characters

            * converting underscores to blanks

            * handling newlines

            * handling comments

            * handling the absence of branch lengths. However, it is optional
              whether the branch length is specified for the root node. If it
              is specified, it is ignored.

        *scalebranchlength* is an optional argument that allows us to rescale
        all branch lengths.  If it has its default value of *None*, then no
        rescaling is done.  If it has another value, then it should be a number
        greater than zero.  In this case, all branch lengths specified by 
        *newick_tree* are multiplied by *scalebranchlength*.

        *tipnames_dict* is an optional argument that allows us to specify
        additional node properties for the tip nodes. If *tipnames_dict* is
        set to a value other than its default of *None*, then it should be
        a dictionary keyed by the node name with a key for each tip node.
        The values are in turn additional dictionaries with key / value pairs.
        For each tip node, the *node.info* dictionary has this key / value
        pair added to it.

        *make_lengths_non_negative* is a Boolean switch specifying that we
        set any branch lengths that are less than zero to be zero
        This is done if the switch is *True*. By default, it is *False*.

        For all tip nodes, the *node.name* attribute is assigned the name string
        given to the node in the Newick tree. Even if the assigned name
        is a number, it is still given as a string.

        Each node of the tree is assigned a unique integer identifier *n*, with
        0 <= *n* < number of nodes in tree.
        """
        assert isinstance(newick_tree, str)
        newick_tree = ''.join(newick_tree.split()) # get rid of any spaces
        assert newick_tree
        assert scalebranchlength == None or (isinstance(scalebranchlength, (int, float)) and scalebranchlength > 0)
        self._newick = newick_tree # string giving Newick format
        # strip away the space and semicolon do some error checking on the Newick tree 
        newick_tree = newick_tree.strip()
        if newick_tree[-1] != ';':
            raise ValueError, "Newick tree does not end in semicolon."
        else: # strip away semicolon and the branch length for the root
            ifinal = len(newick_tree) - 1
            while ifinal > 0 and (newick_tree[ifinal] not in [']', ')']):
                ifinal -= 1
            newick_tree = newick_tree[ : ifinal + 1] 
        if not (newick_tree.count('(') == newick_tree.count(')') > 0):
            raise ValueError, "Invalid number of begin and end parentheses in Newick tree."
        # initialize the tree root
        if newick_tree[-1] == ']': # read info
            j = GetBalancedIndex(newick_tree, len(newick_tree) - 1)
            info = GetInfo(newick_tree[j : ])
            newick_tree = newick_tree[ : j]
        else:
            info = {}
        self._root = Node(root=True, info=info)
        assert newick_tree[0] == '(' and newick_tree[-1] == ')'
        newick_tree = newick_tree[1 : -1]
        RecursivelyAssignDescendents(self._root, newick_tree, scalebranchlength=scalebranchlength, tipnames_dict=tipnames_dict, make_lengths_non_negative=make_lengths_non_negative)
        RecursivelySetNumbers(self._root, 0)
    
    def GetRoot(self):
        """Returns the *Node* object that is the root of the tree."""
        return self._root
    
    def GetNewickTree(self):
        """Returns the Newick string representing the original input tree."""
        return self._newick

    def WriteNewick(self, nodeinfo_keys={}, nodeinfo_res={}, branchinfo_keys={}, branchinfo_res={}):
        """Writes the Newick string corresponding a tree.

        This is NOT the original Newick string used to create the tree (for 
        for that use the *GetNewickTree* method).  Rather, nodes for this 
        tree are labeled according to the arguments *nodeinfo_keys* and 
        *nodeinfo_res* and branches for this tree are formatted according to 
        *branchinfo_keys* and *branchinfo_res*.  
        
        These four arguments are 
        all dictionaries. For any node or branch that has comments
        (specified in the *node.info* dictionary) with a key that
        matches a key in *nodeinfo_keys* or *branchinfo_keys*,
        a comment is printed in the tree. In addition, any comment
        with a key that matches a regular expression that is a key
        in *nodinfo_res* or *branchinfo_res* is printed. The values
        corresponding to the keys can be *None*, in which case the
        method tries to print the correct form of the value (will 
        work for strings, integers, and floats). It can also be
        a function, in which case the method prints the string,
        integer, or float returned by that function which takes
        as input the node or branch information dictionary
        in question. If the function
        returns *None*, then nothing is printed for that value.
        """
        # first define private recursive function
        def _GetString(node):
            """Private recursive function."""
            comments = [] # comments for this node
            for (key, value) in node.info.iteritems():
                assert isinstance(key, str)
                printthis = False
                if key in nodeinfo_keys:
                    printthis = True
                    mapfunc = nodeinfo_keys[key]
                for (m, x) in nodeinfo_res.iteritems():
                    if printthis:
                        break
                    if m.search(key):
                        printthis = True
                        mapfunc = x
                if printthis:
                    if mapfunc != None:
                        value = mapfunc(node.info)
                    if isinstance(value, str):
                        comments.append('%s="%s"' % (key, value))
                    elif isinstance(value, int):
                        comments.append('%s=%d' % (key, value))
                    elif isinstance(value, float):
                        comments.append('%s=%f' % (key, value))
                    else:
                        raise ValueError("'value' has invalid type of %s" % str(type(value)))
            if comments:
                comments = '[&%s]' % (','.join(comments))
            else:
                comments = ''
            if node.name:
                name = "%s" % node.name
            else:
                name = ''
            if node.tip:
                return '%s%s' % (name, comments)
            else:
                # get branch comments
                rightbranchstring = []
                leftbranchstring = []
                for (branchstring, branchinfo) in [(rightbranchstring, node.rightbranchinfo), (leftbranchstring, node.leftbranchinfo)]:
                    for (key, value) in branchinfo.iteritems():
                        assert isinstance(key, str)
                        printthis = False
                        if key in branchinfo_keys:
                            printthis = True
                            mapfunc = branchinfo_keys[key]
                        for (m, x) in branchinfo_res.iteritems():
                            if printthis:
                                break
                            if m.search(key):
                                printthis = True
                                mapfunc = x
                        if printthis:
                            if mapfunc != None:
                                value = mapfunc(branchinfo)
                            if value == None:
                                pass
                            elif isinstance(value, str):
                                branchstring.append('%s="%s"' % (key, value))
                            elif isinstance(value, int):
                                branchstring.append('%s=%d' % (key, value))
                            elif isinstance(value, float):
                                branchstring.append('%s=%f' % (key, value))
                            else:
                                raise ValueError("'value' has invalid type of %s" % str(type(value)))
                if rightbranchstring:
                    rightbranchstring = '[&%s]' % (','.join(rightbranchstring))
                else:
                    rightbranchstring = ''
                if leftbranchstring:
                    leftbranchstring = '[&%s]' % (','.join(leftbranchstring))
                else:
                    leftbranchstring = ''
                rightstring = _GetString(node.rightdescendent)
                leftstring = _GetString(node.leftdescendent)
                return '(%s:%s%f,%s:%s%f)%s%s' % (rightstring, rightbranchstring, node.rightbranch, leftstring, leftbranchstring, node.leftbranch, name, comments)
        # call private recursive function on the tree
        return '%s;' % (_GetString(self._root))



class Node(object):
    """Class for representing nodes of a phylogenetic tree.
    
    Currently only competent for representing nodes on a bifurcating tree.

    The calling arguments for creation of a *Node* object are specified in
    the *__init__* method.

    A node object *n* has the following public attributes:

        * *n.tip* : *True* if the node is a tip, *False* otherwise.

        * *n.root* : *True* if the node is a root node, *False* otherwise.

        * *n.number* : an integer that can be assigned to the node as an identifier.
          It is *None* of no integer has been set.  Typically, this value
          is used to uniquely identify nodes in a tree, and is set during
          construction of the tree.

        * *n.name* : a string giving the name of a node, or *None* if no name exists.

        * *n.info* : a dictionary keyed by arbitrary strings, with values giving
          the corresponding value for that node property.

        * *n.ancestor* : the node that is the ancestor of this node, or *None* if
          the node has no ancestor.  Typically there will be no ancestor
          if *n.root* is *True*.

        * *n.rightdescendent* : node that is the right descendent of this node,
          or *None* if there is no right descendent.  Typically there will be
          no right descendent if *n.tip* is *True*.

        * *n.leftdescendent* : like *n.rightdescendent*, but for the left descendent.

        * *n.ancestorbranch* : the length of the branch to the ancestor, or *None*
          if no such length is set.

        * *n.rightbranch* : the length to the branch to the right descendent, or 
          *None* if no such length is set.

        * *n.leftbranch* : like *n.rightbranch*, but for the left ancestor.

        * *n.ancestorbranchinfo* : a dictionary keyed by arbitrary strings, 
          with values giving the corresponding value for the ancestral branch.

        * *n.rightbranchinfo* : like *n.ancestorbranchinfo*, but for the right branch.

        * *n.leftbranchinfo* : like *n.ancestorbranchinfo*, but for the left branch.

    """
    
    def __init__(self, root=False, tip=False, name=None, ancestor=None, leftdescendent=None, rightdescendent=None, ancestorbranch=None, leftbranch=None, rightbranch=None, number=None, info={}, ancestorbranchinfo={}, rightbranchinfo={}, leftbranchinfo={}):
        """Creates a *Node* object.

        There are no required arguments.  There are the following optional arguments:

        *root* should be *True* iff this is the root of the tree; *False* by default.

        *tip* should be *True* iff this is a tip of the tree; *False* by default.

        *name* can be a string giving the name of the node, if such a name exists.
        By default is *None*.

        *ancestor* can be another *Node* object which is the ancestor of this node.
        By default is *None*, and must be *None* if *root* is *True*.

        *leftdescendent* can be another *Node* object which is the left descendent of 
        this node.  By default is *None*, and must be *None* if *tip* is *True*.

        *rightdescendent* is like *leftdescendent*, but is the right descendent.

        *ancestorbranch* can be a number giving the length of the branch back to the
        ancestor of this node.  By default is *None*, and must be *None* if 
        *ancestor* is *None*.

        *leftbranch* can be a number giving the length of the left branch of the tree.
        By default is *None', and must be *None* if *leftdescendent* is *None*.

        *rightbranch* is like *leftbranch*, but for the right branch.

        *number* is an integer that can be assigned to the node as an identifier.  It
        is *None* by default. 

        *info* is a dictionary that can specify arbitrary information for the node.

        *ancestorbranchinfo* is a dictionary that can specify arbitrary information
        for the ancestor branch.

        *rightbranchinfo* is a dictionary that can specify arbitrary information
        for the right branch.

        *leftbranchinfo* is a dictionary that can specify arbitrary information
        for the left branch.
        """
        assert isinstance(root, bool)
        self.root = root
        assert isinstance(tip, bool)
        self.tip = tip
        assert name == None or isinstance(name, str)
        self.name = name
        assert ancestor == None or (not self.root and isinstance(ancestor, Node))
        self.ancestor = ancestor
        assert leftdescendent == None or (not self.tip and isinstance(leftdescendent, Node))
        assert rightdescendent == None or (not self.tip and isinstance(rightdescendent, Node))
        self.leftdescendent = leftdescendent
        self.rightdescendent = rightdescendent
        assert ancestorbranch == None or isinstance(ancestorbranch, (int, float))
        assert rightbranch == None or isinstance(rightbranch, (int, float))
        assert leftbranch == None or isinstance(leftbranch, (int, float))
        assert not (self.root and ancestorbranch != None)
        assert not (self.tip and rightbranch != None)
        assert not (self.tip and leftbranch != None)
        self.ancestorbranch = ancestorbranch
        self.rightbranch = rightbranch
        self.leftbranch = leftbranch
        assert number == None or isinstance(number, int)
        self.number = number
        assert isinstance(info, dict)
        self.info = info
        assert isinstance(ancestorbranchinfo, dict)
        self.ancestorbranchinfo = ancestorbranchinfo
        assert isinstance(rightbranchinfo, dict)
        self.rightbranchinfo = rightbranchinfo
        assert isinstance(leftbranchinfo, dict)
        self.leftbranchinfo = leftbranchinfo


# test with doctest
if __name__ == '__main__':
    import doctest
    doctest.testmod()
