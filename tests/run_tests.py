"""Runs tests on mutpath package.

Runs doctest on every *.py module file in the package main directory.

Then runs every test in the current directory, which is every file with
a name conforming to test_*.py.

Written by Jesse Bloom.
"""


import os
import sys
import glob
import doctest
import unittest


def main():
    """Runs the tests.
    """

    failurestrings = [] # list of failure strings
    packagename = 'mutpath' # name of package to be tested

    # build_path is the directory with the current build
    build_path = glob.glob('../build/lib*/')
    if (len(build_path) > 1): 
        sys.stderr.write("\nFound more than one build directory." +\
                " Try deleting the build directory and building again.\n")
        raise IOError
    if not build_path:
        sys.stderr.write("\nFailed to find a build directory in ../build" +\
                "\nDid you run 'python setup.py build' before testing?\n")
        raise IOError
    build_path = build_path[0]
    if not os.access(build_path, os.F_OK):
        sys.stderr.write("\nCurrent user does not have access to %s." %\
                build_path + " Try rebuilding as current user.\n")
        raise IOError
    sys.path.insert(1, build_path)
    if not os.path.isdir("%s/%s" % (build_path, packagename)):
        sys.stderr.write("\nCannot find directory with packagename %s" %\
                packagename + ' in build_path %s\n' % build_path) 
        raise IOError

    # Modules to be tested by doctest. 
    # Test everything fitting build_path/*/*.py
    docstring_modnames = ["%s.%s" % (packagename, 
            os.path.splitext(os.path.basename(f))[0]) for f in
            glob.glob('%s/%s/*.py' % (build_path, packagename)) if
            ('__init__' not in f)]
    for modname in docstring_modnames:
        sys.stderr.write("\nTesting %s with doctest... " % (modname))
        module = __import__(modname, None, None, modname.split('.'))
        suite = doctest.DocTestSuite(module)
        del module
        result = unittest.TestResult()
        suite.run(result)
        if result.wasSuccessful():
            sys.stderr.write("all %d tests were successful.\n" %\
                    result.testsRun)
        else:
            sys.stderr.write("test FAILED!\n")
            for (testcase, failstring) in result.failures:
                failurestrings.append(failstring)

    # Test all unittests in test directory with name test_*.py
    for test in glob.glob('test_*.py'):
        test = os.path.splitext(test)[0]
        sys.stderr.write("\nRunning %s... " % test)
        suite = unittest.TestLoader().loadTestsFromName(test)
        result = unittest.TestResult()
        suite.run(result)
        if result.wasSuccessful():
            sys.stderr.write('\nAll tests were successful.\n')
        else:
            sys.stderr.write("\nTest FAILED!\n")
            for (testcase, failstring) in result.failures + result.errors:
                failurestrings.append(failstring)

    # print summary of failures
    if not failurestrings:
        sys.stderr.write("\nTesting complete. All tests were passed successfully.\n")
    else:
        sys.stderr.write("\nTesting complete. Failure on %d the tests."\
                % len(failurestrings))
        sys.stderr.write("\nHere are all of the failures:\n")
        for failstring in failurestrings:
            sys.stderr.write('\n*****************')
            sys.stderr.write("\n%s\n" % failstring)
            sys.stderr.write('\n*****************\n')



if __name__ == '__main__':
    main() # run the program
