"""Module for performing statistics operations operations.


List of functions
--------------------

*MedianCredibleInterval* : Returns posterior median and centered credible interval.


Details of functions
----------------------
Details of functions are listed below in their individual docstrings.
"""


import math
import copy


def MedianCredibleInterval(data, interval):
    """Returns posterior median and centered Bayesian credible interval.
   
    *data* is a list of numbers (must contain at least one number). These
    are assumed to represent numbers sampled from a posterior distribution.

    *interval* is a number with *0 < interval <= 1*.

    This function returns the following tuple *(median, mininterval, maxinterval)*.
    In this tuple:

        * *median* is the median value in *data*. This represents the posterior
          median. 

        * *mininterval* and *maxinterval* are numbers with 
          *maxinterval >= mininterval* such that a fraction *interval* of
          that data in *data* are >= *minterval* and <= *maxinterval*.
          Furthermore, this interval is centered such that the fraction of
          the entries in *data* that are < *minterval* is the same as
          the fraction of entries that are > *maxinterval*. In other words,
          *minterval* and *maxinterval* specify the centered Bayesian credible
          interval.

    >>> data = [0.1, 2, 1, 1.1, 2, 5, -2.1, 7.9, 2, 1.3, 8.1]
    >>> (median, mininterval, maxinterval) = MedianCredibleInterval(data, 0.8)
    >>> median == 2
    True
    >>> mininterval == 0.1
    True
    >>> maxinterval == 7.9
    True
    """
    if not data:
        raise ValueError("No entries")
    if not 0 < interval <= 1:
        raise ValueError("Invalid interval")
    n = len(data)
    datacopy = copy.copy(data)
    datacopy.sort()
    if not (n % 2):
        # even length
        median = (datacopy[n // 2] + datacopy[n // 2 - 1]) / 2.0
    else:
        # odd length
        median = datacopy[n // 2]
    halfoutinterval = (1.0 - interval) / 2.0
    mininterval = datacopy[int(math.floor(halfoutinterval * n))]
    maxinterval = datacopy[int(math.ceil((1.0 - halfoutinterval) * n)) - 1]
    return (median, mininterval, maxinterval)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
