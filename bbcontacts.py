#!/usr/bin/env python
""" bbcontacts is a program to predict residue-level beta-beta contacts
starting from a matrix of predicted residue-residue couplings and a secondary structure assignment or prediction.
"""

import os

import hmm.iohmm
import hmm.paramshmm
import hmm.viterbi
import hmm.postprocesshmm


def main():
    # Keep track of how long this run will take
    import time
    starttime = time.time()

    # Parse options and arguments
    import optparse
    parser = optparse.OptionParser(usage="%prog [options] couplingmatrix diversityvalue outputprefix (-d DSSPsecstructfile | -p PSIPREDsecstructfile)")

    parser.add_option("-d", "--dssp-file", dest="dsspfile", help="Use DSSP file DSSPFILE containing secondary structure assignment")
    parser.add_option("-p", "--psipred-file", dest="psipredfile", help="Use PSIPRED file PSIPREDFILE containing secondary structure prediction")
    parser.add_option("-c", "--config-file", dest="configfile", help="Use config file CONFIGFILE")
    parser.add_option("-s", "--smoothingsize", dest="smoothingsize", type=int, default=10, help="Perform smoothing of the coupling matrix before decoding, over an area extending by SMOOTHINGSIZE in each direction [default=%default, use 0 for no smoothing]")
    parser.add_option("-l", "--long-predictions", "--no-shortening-mode", action="store_true", dest="noshorteningmode", default=False, help="Turn off (slow) prediction-shortening mode (this mode is on by default but will only apply when long predictions occur)")
    parser.add_option("-n", "--pdb-name", dest="pdbname", help="Provide a PDB identifier (when also using -e, this will be the PDB name to look for in EVALUATIONFILE)")
    parser.add_option("-e", "--evaluation-file", dest="evaluationfile", help="Provide a file containing the true contacts (BetaSheet916.dat, BetaSheet1452.dat or same format) for evaluation")

    options, args = parser.parse_args()
    
    # Check that the options and arguments behave as expected
    if len(args) != 3:
        parser.error("Need three positional arguments: coupling matrix, diversity value for that matrix (sqrt(N)/L) and prefix for the output files")
    if (options.dsspfile and options.psipredfile):
        parser.error("Options -d and -p are mutually exclusive")
    if (not options.dsspfile and not options.psipredfile):
        parser.error("One secondary structure file should be provided: DSSP file (-d) or PSIPRED file (-p)")
    if options.evaluationfile and not options.pdbname:
        parser.error("For evaluation, the PDB name should also be provided (-n)")

    couplingmatrix, diversityvalue, outputprefix = args
    diversityvalue = float(diversityvalue)

    if options.evaluationfile and not os.path.exists(options.evaluationfile):
        print "\nWARNING: provided evaluation file %s does not exist, evaluation will not be performed"%options.evaluationfile
        options.evaluationfile = False

    progdir = os.path.dirname(os.path.realpath(__file__)) # current directory
    # Read config file
    from ConfigParser import SafeConfigParser
    config = SafeConfigParser()
    config.optionxform = str # case-sensitive options in config file
    if options.configfile:
        config.read(options.configfile)
    else:
        config.read(os.path.join(progdir, "bbcontacts.conf"))
    
    ## Files containing HMM parameters
    trprobfile = config.get("HMM parameter files", "trprobfile")
    prioroffsetDSSPfile = config.get("HMM parameter files", "prioroffsetDSSPfile")
    prioroffsetPSIPREDfile = config.get("HMM parameter files", "prioroffsetPSIPREDfile")
    condsecstructprobDSSPfile = config.get("HMM parameter files", "condsecstructprobDSSPfile")
    condsecstructprobPSIPREDfile = config.get("HMM parameter files", "condsecstructprobPSIPREDfile")
    singlesecstructprobDSSPfile = config.get("HMM parameter files", "singlesecstructprobDSSPfile")
    singlesecstructprobPSIPREDfile = config.get("HMM parameter files", "singlesecstructprobPSIPREDfile")
    fitparamsfile = config.get("HMM parameter files", "fitparamsfile")

    ## Other constant parameter values
    # Lowest and highest diversity values that the program can manage, given the fits
    lowestdiversity = float(config.get("Diversity parameters", "lowestdiversity"))
    highestdiversity = float(config.get("Diversity parameters", "highestdiversity"))
    # Thresholds for stopping the Viterbi algorithm
    viterbiparams = {}
    section = "Viterbi parameters"
    for option in config.options(section):
        if option in ["tspdbsize", "xtspdbsize"]:
            viterbiparams[option] = int(config.get(section, option))
        else:
            viterbiparams[option] = float(config.get(section, option))
    # Prediction-shortening mode (PSM) parameters
    PSMparams = {}
    section = "PSM parameters"
    for option in config.options(section):
        if option=="PSMtrprobstepsize":
            PSMparams[option] = float(config.get(section, option))
        else:
            PSMparams[option] = int(config.get(section, option))
    # Masking parameters
    maskarounddiagparallel = int(config.get("Masking parameters", "maskarounddiagparallel"))
    maskarounddiagantiparallel = int(config.get("Masking parameters", "maskarounddiagantiparallel"))
    maskaroundcontact = int(config.get("Masking parameters", "maskaroundcontact"))

    # Identifier for this run
    if options.pdbname:
        identifier = options.pdbname
    else:
        # Generate an identifier from the coupling matrix filename
        pdb = os.path.basename(couplingmatrix).split(".")[0]
        if len(pdb)==9:
            identifier = pdb[0:4]+pdb[5]+pdb[7:9]
        else:
            identifier = pdb

    print "\nStarting run %s %.3f"%(identifier,diversityvalue)

    # Process the secondary structure input file
    dsspMasking = False
    if options.dsspfile:
        # Should we apply DSSP masking?
        dsspMasking = True
        print "\nWARNING: DSSP masking is on!!"
        secstructseq, secstructdic = hmm.iohmm.processSecStruct(options.dsspfile, outputprefix, identifier, diversityvalue, dsspMasking, options.evaluationfile)
        outputprefix += ".DSSP"
    else:
        dsspMasking = False
        secstructseq, secstructdic = hmm.iohmm.processSecStruct(options.psipredfile, outputprefix, identifier, diversityvalue, dsspMasking, options.evaluationfile)

    # Retrieve couplings
    predcouplings, nbres = hmm.iohmm.retrieveCouplings(couplingmatrix, outputprefix, identifier, diversityvalue, secstructseq, options.evaluationfile)
    if options.smoothingsize and options.smoothingsize != 0:
        # Smooth the coupling matrix
        predcouplings = hmm.paramshmm.smoothMatrix(predcouplings, nbres, options.smoothingsize)

    # Retrieve parameters for the secondary-structure-based part of the emission probabilities
    if dsspMasking:
        secprobdic = hmm.paramshmm.getSecondaryStructureProbabilities(condsecstructprobDSSPfile)
        secprior = hmm.paramshmm.getSecondaryStructurePrior(singlesecstructprobDSSPfile)
    else:
        secprobdic = hmm.paramshmm.getSecondaryStructureProbabilities(condsecstructprobPSIPREDfile)
        secprior = hmm.paramshmm.getSecondaryStructurePrior(singlesecstructprobPSIPREDfile)

    # Retrieve parameters for the coupling-based part of the emission probabilities
    if diversityvalue > highestdiversity:
        print "\nWARNING - Diversity too high, setting diversity to %.3f"%highestdiversity
        diversityvalue = highestdiversity
    params, lowDiversity = hmm.paramshmm.getTransformedGammaParameters(identifier, diversityvalue, fitparamsfile)
    if diversityvalue < lowestdiversity:
        lowDiversity = True
    if lowDiversity:
        print "\nWARNING: there are too few sequences in the multiple sequence alignment, so only the secondary structure will count to predict beta contacts (residue couplings will be ignored) for %s %.2f"%(identifier,diversityvalue)

    # Retrieve transition probabilities
    trprob = hmm.paramshmm.getTransitionProbabilities(trprobfile)

    # Retrieve prior offset (depending on distance to diagonal)
    if dsspMasking:
        prioroffset = hmm.paramshmm.getPriorOffset(prioroffsetDSSPfile,nbres)
    else:
        prioroffset = hmm.paramshmm.getPriorOffset(prioroffsetPSIPREDfile,nbres)

    # Set thresholds for stopping the Viterbi algorithm
    viterbirecalclower, viterbirecalclowerPSM = hmm.paramshmm.setViterbiThresholds(identifier, nbres, dsspMasking, viterbiparams)
    
    # Calculate all emission probabilites
    print "\nCalculating emission probabilities"
    emissions = hmm.paramshmm.calculateAllEmissions(nbres, predcouplings, secstructseq, secprobdic, params, lowDiversity)

    # Run the Viterbi algorithm
    print "\nViterbi paths"
    allpaths, recalculateCount = hmm.viterbi.runViterbi(nbres, trprob, prioroffset, secprior, secprobdic, emissions, viterbirecalclower, viterbirecalclowerPSM, PSMparams, dsspMasking, secstructseq, secstructdic, maskarounddiagparallel, maskarounddiagantiparallel, maskaroundcontact, options.noshorteningmode)

    # Post-process results and write output
    print "\nPost-processing results and writing output..."
    hmm.postprocesshmm.postProcessPaths(allpaths,outputprefix,viterbirecalclower,identifier,diversityvalue, maskaroundcontact, secstructdic, options.evaluationfile)

    endtime = time.time()

    print "End of run: %s %5.3f %6i\tPSM iterations: %3i\tLength of run: %16.2f s"%(identifier, diversityvalue, nbres, recalculateCount-1, endtime-starttime)

if __name__ == "__main__":
    main()
