"""Module for functions related to reading / writing BEAST files.

This module is designed to interface with the BEAST phylogenetics
package.

Written by Jesse Bloom.
"""


import os



def WriteXML(xmlfile, seqs, datatype, chainlength, screenloginterval,\
        fileloginterval, usemarkovjumps, sitemodel):
    """Writes an XML file that can serve as input to BEAST.

    Dates are assigned to sequences based on the values specified
    in the names in seqs.

    The file logs are given the same base prefix as that for *xmlfile* itself,
    but the extensions are .log and .trees.

    *xmlfile* a string giving the name of the XML file that
    we create. It is overwritten if it already exists.

    *seqs* a list of sequences in the form (header, seq). 
    These sequences are assumed to have dates assigned,
    and are specified following a terminal underscore,
    such as::

        A/Nanchang/933/1995_1995.00

    *datatype* the data type for the sequences, might be
    'nucleotide' or 'amino acid'.  

    *chainlength* integer giving the chain length for the MCMC.

    *screenloginterval* integer specifying the frequency with which
    we log MCMC information to the screen.

    *fileloginterval* integer specifying the frequency with which
    we log MCMC information to the file logs.

    *usemarkovjumps* Boolean switch. If True, we use the Markov
    Jumps feature to write mutations / sequences for the .trees
    file. This allows reconstruction of mutational paths, but
    also increases the .trees file size dramatically.

    *sitemodel* the site substitution model. Currently allows
    'JTT' for proteins and 'HKY' for nucleic acids.

    Various other options are currently hardcoded in the method.
    Options could be set to allow them to be changed. Currently:

    * uses a coalescent-based prior on the trees with constant
      population size with Jeffrey's prior.

    * uses a strict clock with uniform rates across branches, gamma 
      prior on rate.

    """
    if sitemodel == 'JTT':
        if datatype != 'amino acid':
            raise ValueError("sitemodel of JTT requires datatype of amino acid")
    elif sitemodel == 'HKY':
        if datatype != 'nucleotide':
            raise ValueError("sitemodel of HKY requires datatype of nucleotide")
    else:
        raise ValueError("Invalid sitemodel %s" % sitemodel)
    basename = os.path.splitext(xmlfile)[0]
    f = open(xmlfile, 'w')
    #
    # write header n
    f.write('<?xml version="1.0" standalone="yes"?>')
    #
    # begin the beast block
    f.write('\n\n<beast>\n\n')
    #
    # write the taxa block with dates
    f.write('\t<taxa id="taxa">\n')
    for (head, seq) in seqs:
        try:
            date = float(head.split('_')[-1])
        except ValueError:
            raise ValueError("problem getting date after underscore in %s" % head)
        f.write('\t\t<taxon id="%s">\n\t\t\t<date value="%.2f" direction="forwards" units="years"/>\n\t\t</taxon>\n' % (head, date))
    f.write('</taxa>\n\n')
    #
    # write alignment block
    f.write('\t<alignment id="alignment" dataType="%s">\n' % datatype)
    for (head, seq) in seqs:
        f.write('\t\t<sequence>\n\t\t\t<taxon idref="%s"/>\n\t\t\t%s\n\t\t</sequence>\n' % (head, seq))
    f.write('\t</alignment>\n')
    #
    # write the patterns block, do not restrict to non-unique sites
    # if using Markov Jumps
    f.write(\
        '\n\t<patterns id="patterns" from="1" unique="%s">\n' % {False:'true', True:'false'}[usemarkovjumps] +\
        '\t\t<alignment idref="alignment"/>\n' +\
        '\t</patterns>\n')
    # 
    # write blocks related to the coalescent tree
    f.write('\n\t<constantSize id="constant" units="years">\n\t\t<populationSize>\n\t\t\t<parameter id="constant.popSize" value="420.0" lower="0.0" upper="Infinity"/>\n\t\t</populationSize>\n\t</constantSize>\n')
    f.write('\n\t<coalescentTree id="startingTree">\n\t\t<taxa idref="taxa"/>\n\t\t<constantSize idref="constant"/>\n\t</coalescentTree>\n')
    f.write('\n\t<treeModel id="treeModel">\n\t\t<coalescentTree idref="startingTree"/>\n\t\t<rootHeight>\n\t\t\t<parameter id="treeModel.rootHeight"/>\n\t\t</rootHeight>\n\t\t<nodeHeights internalNodes="true">\n\t\t\t<parameter id="treeModel.internalNodeHeights"/>\n\t\t</nodeHeights>\n\t\t<nodeHeights internalNode="true" rootNode="true">\n\t\t\t<parameter id="treeModel.allInternalNodeHeights"/>\n\t\t</nodeHeights>\n\t</treeModel>\n')
    f.write('\n<coalescentLikelihood id="coalescent">\n\t\t<model>\n\t\t\t<constantSize idref="constant"/>\n\t\t</model>\n\t\t<populationTree>\n\t\t\t<treeModel idref="treeModel"/>\n\t\t</populationTree>\n\t</coalescentLikelihood>\n')
    #
    # write strict clock block
    f.write(\
        '\n\t<strictClockBranchRates id="branchRates">\n' +\
        '\t\t<rate>\n\t\t\t<parameter id="clock.rate" value="7.0e-5" lower="0.0" upper="Infinity"/>\n' +\
        '\t\t</rate>\n' +\
        '\t</strictClockBranchRates>\n')
    #
    # write the JTT or HKY site model block
    if sitemodel == 'JTT':
        f.write('\n\t<aminoAcidModel id="aa" type="JTT"/>\n')
        f.write(\
            '\n\t<siteModel id="siteModel">\n' +\
            '\t\t<substitutionModel>\n' +\
            '\t\t\t<aminoAcidModel idref="aa"/>\n' +\
            '\t\t</substitutionModel>\n' +\
            '\t</siteModel>\n')
    elif sitemodel == 'HKY':
        f.write(\
             '\n\t<HKYModel id="hky">\n' +\
             '\t\t<frequencies>\n' +\
             '\t\t\t<frequencyModel dataType="nucleotide">\n' +\
             '\t\t\t\t<frequencies>\n' +\
             '\t\t\t\t\t<parameter id="frequencies" value="0.25 0.25 0.25 0.25"/>\n' +\
             '\t\t\t\t</frequencies>\n' +\
             '\t\t\t</frequencyModel>\n' +\
             '\t\t</frequencies>\n' +\
             '\t\t<kappa>\n' +\
             '\t\t\t<parameter id="kappa" value="2.0" lower="0.0" upper="Infinity"/>\n' +\
             '\t\t</kappa>\n' +\
             '\n\t</HKYModel>\n')
        f.write(\
            '\n\t<siteModel id="siteModel">\n' +\
            '\n\t\t<substitutionModel>\n' +\
            '\n\t\t\t<HKYModel idref="hky"/>\n' +\
            '\n\t\t</substitutionModel>\n' +\
            '\n\t</siteModel>\n')
    else:
        raise ValueError("Invalid sitemodel of %s" % sitemodel)
    #
    # write either the standard or Markov Jumps tree-likelihood block
    if usemarkovjumps:
        f.write(\
            '\n\t<markovJumpsTreeLikelihood id="treeLikelihood" useUniformization="true" saveCompleteHistory="true" logCompleteHistory="true" compactHistory="true">\n' +\
            '\t\t<patterns idref="patterns"/>\n' +\
            '\t\t<treeModel idref="treeModel"/>\n' +\
            '\t\t<siteModel idref="siteModel"/>\n')
        if sitemodel == 'JTT':
            f.write('\t\t<aminoAcidModel idref="aa"/>\n')
        elif sitemodel == 'HKY':
            f.write('\t\t<HKYModel idref="hky"/>\n')
        else:
            raise ValueError("Invalid sitemodel of %s" % sitemodel)
        f.write(\
            '\t\t<strictClockBranchRates idref="branchRates"/>\n' +\
            '\t</markovJumpsTreeLikelihood>\n')
    else:
        f.write(\
            '\n\t<treeLikelihood id="treeLikelihood" useAmbiguities="false">\n' +\
            '\t\t<patterns idref="patterns"/>\n' +\
            '\t\t<treeModel idref="treeModel"/>\n' +\
            '\t\t<siteModel idref="siteModel"/>\n' +\
            '\t\t<strictClockBranchRates idref="branchRates"/>\n' +\
            '\t</treeLikelihood>\n')

    #
    # write operators block
    f.write(\
        '\n\t<operators id="operators" optimizationSchedule="log">\n')
    if sitemodel == 'HKY':
        f.write(\
            '\t\t<scaleOperator scaleFactor="0.75" weight="0.1">\n' +\
            '\t\t\t<parameter idref="kappa"/>\n' +\
            '\t\t</scaleOperator>\n' +\
            '\t\t<deltaExchange delta="0.01" weight="0.1">\n' +\
            '\t\t\t<parameter idref="frequencies"/>\n' +\
            '\t\t</deltaExchange>\n')
    elif sitemodel == 'JTT':
        pass
    else:
        raise ValueError("Invalid sitemodel of %s" % sitemodel)
    f.write(\
        '\t\t<scaleOperator scaleFactor="0.75" weight="3">\n' +\
        '\t\t\t<parameter idref="clock.rate"/>\n' +\
        '\t\t</scaleOperator>\n' +\
        '\t\t<subtreeSlide size="43.0" gaussian="true" weight="15">\n' +\
        '\t\t\t<treeModel idref="treeModel"/>\n' +\
        '\t\t</subtreeSlide>\n' +\
        '\t\t<narrowExchange weight="15">\n' +\
        '\t\t\t<treeModel idref="treeModel"/>\n' +\
        '\t\t</narrowExchange>\n' +\
        '\t\t<wideExchange weight="3">\n' +\
        '\t\t\t<treeModel idref="treeModel"/>\n' +\
        '\t\t</wideExchange>\n' +\
        '\t\t<wilsonBalding weight="3">\n' +\
        '\t\t\t<treeModel idref="treeModel"/>\n' +\
        '\t\t</wilsonBalding>\n' +\
        '\t\t<scaleOperator scaleFactor="0.75" weight="3">\n' +\
        '\t\t\t<parameter idref="treeModel.rootHeight"/>\n' +\
        '\t\t</scaleOperator>\n' +\
        '\t\t<uniformOperator weight="30">\n' +\
        '\t\t\t<parameter idref="treeModel.internalNodeHeights"/>\n' +\
        '\t\t</uniformOperator>\n' +\
        '\t\t<scaleOperator scaleFactor="0.75" weight="3">\n' +\
        '\t\t\t<parameter idref="constant.popSize"/>\n' +\
        '\t\t</scaleOperator>\n' +\
        '\t\t<upDownOperator scaleFactor="0.75" weight="3">\n' +\
        '\t\t\t<up>\n' +\
        '\t\t\t\t<parameter idref="clock.rate"/>\n' +\
        '\t\t\t</up>\n' +\
        '\t\t\t<down>\n' +\
        '\t\t\t\t<parameter idref="treeModel.allInternalNodeHeights"/>\n' +\
        '\t\t\t</down>\n' +\
        '\t\t</upDownOperator>\n' +\
        '\t</operators>\n')
    # 
    # write the mcmc block
    f.write(\
        '\n\t<mcmc id="mcmc" chainLength="%d" autoOptimize="true" autoOptimizeDelay="10000">' % chainlength +\
        '\t\t<posterior id="posterior">\n' +\
        '\t\t\t<prior id="prior">\n')
    if sitemodel == 'HKY':
        f.write(\
            '\t\t\t\t<logNormalPrior mean="1.0" stdev="1.25" offset="0.0" meanInRealSpace="false">\n' +\
            '\t\t\t\t\t<parameter idref="kappa"/>\n' +\
            '\t\t\t\t</logNormalPrior>\n')
    elif sitemodel == 'JTT':
        pass
    else:
        raise ValueError("invalid sitemodel %s" % sitemodel)
    f.write(\
        '\t\t\t\t<oneOnXPrior>\n' +\
        '\t\t\t\t\t<parameter idref="constant.popSize"/>\n' +\
        '\t\t\t\t</oneOnXPrior>\n' +\
        '\t\t\t\t<gammaPrior shape="1.0" scale="1.0" offset="0.0">\n' +\
        '\t\t\t\t\t<parameter idref="clock.rate"/>\n' +\
        '\t\t\t\t</gammaPrior>\n' +\
        '\t\t\t\t<coalescentLikelihood idref="coalescent"/>\n' +\
        '\t\t\t</prior>\n' +\
        '\t\t\t<likelihood id="likelihood">\n' +\
        '\t\t\t\t<markovJumpsTreeLikelihood idref="treeLikelihood"/>\n' +\
        '\t\t\t</likelihood>\n' +\
        '\t\t</posterior>\n' +\
        '\n\t\t<operators idref="operators"/>\n' +\
        '\n\t\t<log id="screenLog" logEvery="%d">\n' % screenloginterval +\
        '\t\t\t<column label="Posterior" dp="4" width="12">\n' +\
        '\t\t\t\t<likelihood idref="posterior"/>\n' +\
        '\t\t\t</column>\n' +\
        '\t\t\t<column label="Prior" dp="4" width="12">\n' +\
        '\t\t\t\t<likelihood idref="prior"/>\n' +\
        '\t\t\t</column>\n' +\
        '\t\t\t<column label="Likelihood" dp="4" width="12">\n' +\
        '\t\t\t\t<likelihood idref="likelihood"/>\n' +\
        '\t\t\t</column>\n' +\
        '\t\t\t<column label="rootHeight" sf="6" width="12">\n' +\
        '\t\t\t\t<parameter idref="treeModel.rootHeight"/>\n' +\
        '\t\t\t</column>\n' +\
        '\t\t\t<column label="clock.rate" sf="6" width="12">\n' +\
        '\t\t\t\t<parameter idref="clock.rate"/>\n' +\
        '\t\t\t</column>\n' +\
        '\t\t</log>\n' +\
        '\n\t\t<log id="fileLog" logEvery="%d" fileName="%s.log" overwrite="false">\n' % (fileloginterval, basename) +\
        '\t\t\t<posterior idref="posterior"/>\n' +\
        '\t\t\t<prior idref="prior"/>\n' +\
        '\t\t\t<likelihood idref="likelihood"/>\n' +\
        '\t\t\t<parameter idref="treeModel.rootHeight"/>\n' +\
        '\t\t\t<parameter idref="constant.popSize"/>\n')
    if sitemodel == 'HKY':
        f.write(\
            '\t\t\t<parameter idref="kappa"/>\n' +\
            '\t\t\t<parameter idref="frequencies"/>\n')
    elif sitemodel == 'JTT':
        pass
    else:
        raise ValueError("invalid sitemodel %s" % sitemodel)
    f.write(\
        '\t\t\t<parameter idref="clock.rate"/>\n' +\
        '\t\t\t<coalescentLikelihood idref="coalescent"/>\n' +\
        '\t\t</log>\n' +\
        '\n\t\t<logTree id="treeFileLog" logEvery="%d" nexusFormat="true" fileName="%s.trees" sortTranslationTable="true">\n' % (fileloginterval, basename) +\
        '\t\t\t<treeModel idref="treeModel"/>\n' +\
        '\t\t\t<posterior idref="posterior"/>\n' +\
        '\t\t\t<markovJumpsTreeLikelihood idref="treeLikelihood"/>\n' +\
        '\t\t</logTree>\n' +\
        '\t</mcmc>\n')
    #
    # report block
    f.write(\
        '\n\t<report>\n' +\
        '\t\t<property name="timer">\n' +\
        '\t\t\t<object idref="mcmc"/>\n' +\
        '\t\t</property>\n' +\
        '\t</report>\n')
    # 
    # end the beast block
    f.write('</beast>')

    f.close()
