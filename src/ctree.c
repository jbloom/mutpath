// Implements fast C versions of some of the functions in the tree
// module for the ctree package. Written by Jesse Bloom.
//
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Documentation string for GetBalancedIndex function
static char GetBalancedIndex_docs[] = "Fast C function findes the index of the balancing parentheses character.\n\nExactly mimics the GetBalancedIndex function in mutpath.tree.";

// Implementation of GetBalancedIndex function
static PyObject *GetBalancedIndex(PyObject *self, PyObject *args) {
    // Calling variables: s (string), i (int)
    char *s;
    char si, partner;
    int lookforward;
    long int i, j, nested;
    size_t lens;
    if (! PyArg_ParseTuple(args, "sl", &s, &i)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to GetBalancedIndex.");
        return NULL;
    }
    lens = strlen(s);
    if ((i < 0) || (i >= lens)) {
        PyErr_SetString(PyExc_ValueError, "Invalid index for string.");
        return NULL;
    }
    si = s[i];
    switch(si)
    {
        case '(':
            partner = ')';
            lookforward = 1;
            break;
        case ')':
            partner = '(';
            lookforward = 0;
            break;
        case '[':
            partner = ']';
            lookforward = 1;
            break;
        case ']':
            partner = '[';
            lookforward = 0;
            break;
        case '{':
            partner = '}';
            lookforward = 1;
            break;
        case '}':
            partner = '{';
            lookforward = 0;
            break;
        default:
            PyErr_SetString(PyExc_ValueError, "Character is not a bracket");
            return NULL;
    }
    nested = 0;
    if (lookforward) {
        j = i + 1;
        while (j < lens) {
            if (s[j] == si) {
                nested++;
            }
            else if (s[j] == partner) {
                if (nested > 0) {
                    nested--;
                }
                else {
                    return PyInt_FromLong(j);
                }
            }
            j++;
        }
        PyErr_SetString(PyExc_ValueError, "Failed to balance");
        return NULL;
    }
    else {
        j = i - 1;
        while (j >= 0) {
            if (s[j] == si) {
                nested++;
            }
            else if (s[j] == partner) {
                if (nested > 0) {
                    nested--;
                }
                else {
                    return PyInt_FromLong(j);
                }
            }
            j--;
        }
        PyErr_SetString(PyExc_ValueError, "Failed to balance");
        return NULL;
    }
}

// Module documentation string
static char ctree_docs[] = "Fast implementations of some functions from mutpath.tree.\n\nGetBalancedIndex in this module mimics the same function from tree.";

// The module functions
static PyMethodDef ctree_funcs[] = {
    {"GetBalancedIndex", (PyCFunction) GetBalancedIndex, METH_VARARGS, GetBalancedIndex_docs},
    {NULL}
};

// Initialize the module
void initctree(void) {
    Py_InitModule3("ctree", ctree_funcs, ctree_docs);
}
