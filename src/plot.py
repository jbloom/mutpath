"""Module for performing plotting for ``mutpath`` package.

This module uses ``pylab`` and ``matplotlib`` to make plots. These plots will
fail if ``pylab`` and ``matplotlib`` are not available for importation. Before
running any function in this module, you can run the *PylabAvailable*
function to determine if ``pylab`` and ``matplotlib`` are available. Otherwise,
calling any other function will raise an Exception if thise modules are
not available. The ``pdf`` backend is used for ``matplotlib`` / ``pylab``. This means
that plots must be created as PDF files.

Functions are:

`PylabAvailable`

`CumulativeFractionPlot`

'DatesPlot`

`Base10Formatter`

`SplitLabel`

Written by Jesse Bloom.
"""


import os
import sys
import math


# global variable _pylabavailable indicates if pylab/matplotlib present
try:
    import matplotlib
    matplotlib.use('pdf')
    import pylab
    _pylabavailable = True
except ImportError:
    _pylabavailable = False


def PylabAvailable():
    """Returns True if pylab/matplotlib available, False otherwise.
    
    You should call this function to test for the availability of the
    pylab/matplotlib plotting modules before using other functions in
    this module.
    """
    return _pylabavailable



def DatesPlot(mutdates, plotfile, interval):
    """Plots dates of mutations.

    Uses pylab / matplotlib to plot the dates and credible intervals
    for mutations. Will raise an error *PylabAvailable() == False*.
    The plot is a PDF.    

    * *mutdates* is a list of the mutations, in the form of the tuples
      *(median, mininterval, maxinterval, mut, fractoca, weight)*. Mutations
      are plotted in the order they are listed. In these tuples:

      * *median* : posterior median date

      * *minterval* : minimum of credible interval

      * *maxinterval* : maximum of credible interval

      * *mut* : string giving name of mutation

      * *fractoca* : probability mutation is on path from common ancestor
        to starting sequence

      * *weight* : fraction of paths containing mutation.
    
    * *plotfile* is a string giving the name of the PDF file we create.

    * *interval* is the range of the credible interval. For example, 0.9
      means a 90% credible interval.
    """
    ext = os.path.splitext(plotfile)[1].lower()
    if ext != '.pdf':
        raise ValueError("Extension must be .pdf, but found %s" % ext)
    if not PylabAvailable():
        raise ValueError("pylab / matplotlib not available.")
    if not mutdates:
        raise ValueError("no mutation dates to plot")
    tocalabels = []
    tocamedians = []
    tocaerrlow = []
    tocaerrhigh = []
    tocays = []
    fromcalabels = []
    fromcamedians = []
    fromcaerrlow = []
    fromcaerrhigh = []
    fromcays = []
    y = 0
    for (median, mininterval, maxinterval, mut, fractoca, weight) in mutdates:
        label = "%s" % (mut)
        errlow = median - mininterval
        errhigh = maxinterval - median
        if fractoca > 0.5:
            tocays.append(y)
            tocalabels.append(label)
            tocamedians.append(median)
            tocaerrlow.append(errlow)
            tocaerrhigh.append(errhigh)
        else:
            fromcays.append(y)
            fromcalabels.append(label)
            fromcamedians.append(median)
            fromcaerrlow.append(errlow)
            fromcaerrhigh.append(errhigh)
        y += 1
    (lmargin, rmargin, bmargin, tmargin) = (0.11, 0.05, 0.08, 0.01)
    matplotlib.rc('font', size=10)
    matplotlib.rc('xtick', labelsize=10) 
    matplotlib.rc('ytick', labelsize=10) 
    matplotlib.rc('legend', numpoints=1)
    matplotlib.rc('legend', fontsize=10)
    fig = pylab.figure(figsize=(6, 6))
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 - tmargin - bmargin])
    tocabar = fromcabar = None
    if tocalabels:
        tocabar = pylab.errorbar(tocamedians, tocays, xerr=[tocaerrlow, tocaerrhigh], fmt='sr')
    if fromcalabels:
        fromcabar = pylab.errorbar(fromcamedians, fromcays, xerr=[fromcaerrlow, fromcaerrhigh], fmt='sb')
    ny = len(mutdates)
    pylab.gca().set_ylim((-1, ny))
    pylab.gca().yaxis.set_major_locator(matplotlib.ticker.FixedLocator([y for y in range(ny)]))
    pylab.gca().yaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(tocalabels + fromcalabels))
    pylab.gca().xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
    pylab.xlabel("Date (posterior median and Bayesian %.2f%s credible interval)" % (interval * 100, '%'))
    if tocabar and fromcabar:
        pylab.legend([tocabar[0], fromcabar[0]], ['path to common ancestor', 'path from common ancestor'], loc='lower right')
    elif tocabar:
        pylab.legend([tocabar[0]], ['path to common ancestor'], loc='lower right')
    elif fromcabar:
        pylab.legend([fromcabar[0]], ['path from common ancestor'], loc='lower right')
    pylab.savefig(plotfile)



def CumulativeFractionPlot(datalist, plotfile, title, xlabel):
    """Creates a cumulative fraction plot.

    Takes a list of numeric data. Plots a cumulative fraction
    plot giving the fraction of the data points that are <=
    the indicated value.

    *datalist* is a list of numbers giving the data for which we
    are computing the cumulative fraction plot. Raises an 
    exception if this is an empty list.

    *plotfile* is the name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be '.pdf'.

    *title* is a string placed above the plot as a title. Uses LaTex
    formatting.

    *xlabel* is the label given to the X-axis. Uses LaTex formatting.

    This function uses pylab / matplotlib. It will raise an Exception if
    these modules cannot be imported (if PylabAvailable() is False).
    """
    if len(datalist) < 1:
        raise ValueError("datalist is empty")
    if not _pylabavailable:
        raise ImportError("Could not find pylab or matplotlib")
    if os.path.splitext(plotfile)[1] != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    datalist.sort() # sort from smallest to largest
    (xmin, xmax) = (datalist[0], datalist[-1])
    n = len(datalist)
    cumfracs = []
    cf = 0.0
    for x in datalist:
        cf += 1. / n
        cumfracs.append(cf)
    assert len(datalist) == len(cumfracs)
    assert abs(1.0 - cf) < 1e-7
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=12)
    fig = pylab.figure(figsize=(6, 4))
    (lmargin, rmargin, bmargin, tmargin) = (0.1, 0.01, 0.15, 0.1)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
            bmargin - tmargin])
    pylab.plot(datalist, cumfracs, 'r-')
    pylab.gca().set_ylim([0, 1])
    pylab.gca().set_xlim([xmin, xmax])
    pylab.ylabel('cumulative fraction')
    pylab.xlabel(xlabel)
    pylab.title(title)
    if plotfile:
        pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def Base10Formatter(number, exp_cutoff, exp_decimal_digits, decimal_digits):
    """Converts a number into Latex formatting with scientific notation.

    Takes a number and converts it to a string that can be shown
    in LaTex using math mode. It is converted to scientific notation
    if the criteria specified by exp_cutoff.

    *number* the number to be formatted, should be a float or integer.
    Currently only works for numbers >= 0

    *exp_cutoff* convert to scientific notation if abs(math.log10(number)) >= this.

    *exp_decimal_digits* show this many digits after the decimal if number
    is converted to scientific notation.

    *decimal_digits* show this many digits after the decimal if number
    is NOT converted to scientific notation.
    
    The returned value is the LaTex' string. If the number is zero, the
    returned string is simply '0'.

    >>> Base10Formatter(103, 3, 1, 1)
    '103.0'

    >>> Base10Formatter(103.0, 2, 1, 1)
    '1.0 \\\\times 10^{2}'

    >>> Base10Formatter(103.0, 2, 2, 1)
    '1.03 \\\\times 10^{2}'

    >>> Base10Formatter(2892.3, 3, 1, 1) 
    '2.9 \\\\times 10^{3}'

    >>> Base10Formatter(0.0, 3, 1, 1) 
    '0'

    >>> Base10Formatter(0.012, 2, 1, 1)
    '1.2 \\\\times 10^{-2}'
  
    >>> Base10Formatter(-0.1, 3, 1, 1)
    Traceback (most recent call last):
        ...
    ValueError: number must be >= 0
    """
    if number < 0:
        raise ValueError('number must be >= 0')
    if number == 0:
        return '0'
    exponent = int(math.log10(number))
    if math.log10(number) < exponent and number < 1:
        exponent -= 1
    if abs(exponent) >= exp_cutoff:
        x = number / (10.**exponent)
        formatstr = '%.' + '%d' % exp_decimal_digits + 'f \\times 10^{%d}'
        return formatstr % (x, exponent)
    else:
        formatstr = '%.' + '%d' % decimal_digits + 'f'
        return formatstr % number


def SplitLabel(label, splitlen, splitchar):
    """Splits a string with a return if it exceeds a certain length.

    *label* a string giving the label we might split.

    *splitlen* the maximum length of a label before we attempt to 
    split it.

    *splitchar* the character added when splitting a label.

    If len(*label*) > *splitlen*, we attempt to split the label in the
    middle by adding *splitchar*. The label is split as close to the
    middle as possible while splitting at a space.

    No splitting as label length less than *splitlen*

    >>> SplitLabel('WT virus 1', 10, '\\n')
    'WT virus 1'

    Splitting of this label

    >>> SplitLabel('WT plasmid 1', 10, '\\n')
    'WT\\nplasmid 1'

    Splitting of this label

    >>> SplitLabel('mutated WT plasmid 1', 10, '\\n')
    'mutated WT\\nplasmid 1'

    """
    if len(label) <= splitlen:
        return label
    else:
        j = 0
        imid = len(label) // 2
        index = None
        while 0 <= imid - j <= imid + j < len(label):
            if label[imid - j].isspace():
                return "%s%s%s" % (label[ : imid - j], splitchar, label[imid - j + 1 : ])
            elif label[imid + j].isspace():
                return "%s%s%s" % (label[ : imid + j], splitchar, label[imid + j + 1 : ])
            j += 1
        else:
            return label # no white space to split


if __name__ == '__main__':
    import doctest
    doctest.testmod()
