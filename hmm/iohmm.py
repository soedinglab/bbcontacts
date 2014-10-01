""" This is a collection of IO scripts for the bbcontacts methods """

import os
import sys
import collections
import re

def processSecStruct(secstructfile, outputprefix, identifier, diversityvalue, dsspMasking, evaluationfile):
    """ Process the secondary structure input file to retrieve a sequence of secstruct states"""
    ssseq = ""
    if not os.path.exists(secstructfile):
        writeOutputBeforeExiting(outputprefix, identifier, diversityvalue, evaluationfile)
        sys.exit("The secondary structure file %s does not exist!"%secstructfile)
    with open(secstructfile, "r") as f:
        for l in f:
            if l[0]==">":
                continue
            ssseq += l.strip()
    # Check that the secondary structure contains only H, E and C
    pattern = re.compile(r'[^HEC]+')
    if bool(pattern.search(ssseq)):
        writeOutputBeforeExiting(outputprefix, identifier, diversityvalue, evaluationfile)
        sys.exit("The secondary structure file %s should contain only H, E and C!"%secstructfile)
    # Add a final C state that will be used for bordering cells of the coupling matrix (for the conditional probability term in the termination of the Viterbi)
    ssseq += "C"
    # The ssdic is used when the input file is a DSSP assignment file to record the positions that are part of DSSP strands
    ssdic = {}
    ssdic["all"] = []
    ssdic["allprev"] = [] # used in Viterbi initialization (start state)
    ssdic["allprevanti"] = [] # used in Viterbi initialization (start state)
    if dsspMasking:
        strandnumber = 0
        for (res, ss) in enumerate(ssseq):
            if ss=="E":
                ssdic[res+1] = strandnumber # for all residues part of the same strand, this number will be the same (used in Viterbi initialization and when filtering the output)
                ssdic["all"].append(int(res)+1)
                ssdic["allprev"].append(int(res)+1-1)
                ssdic["allprevanti"].append(int(res)+1+1)
            else:
                strandnumber += 1
    return ssseq, ssdic


def retrieveCouplings(couplingmatrix, outputprefix, identifier, diversityvalue, secstructseq, evaluationfile):
    """ Retrieve couplings from the input couplingmatrix file """
    pred = {}
    if not os.path.exists(couplingmatrix):
        writeOutputBeforeExiting(outputprefix, identifier, diversityvalue, evaluationfile)
        sys.exit("The coupling matrix file %s does not exist!"%couplingmatrix)

    ff = open(couplingmatrix).readlines()
    for ite,l in enumerate(ff):
        # start indices from 1 and increment by 1
        pred[ite+1] = {}
        for ite2,m in enumerate(l.strip().split("\t")):
            pred[ite+1][ite2+1] = float(m) # matrix of prediction scores

    nbres = max(pred) # number of residues in this protein
    if nbres != len(secstructseq)-1:
        writeOutputBeforeExiting(outputprefix, identifier, diversityvalue, evaluationfile)
        sys.exit("Problem with %s: The size of the coupling matrix does not match the length of the secondary structure sequence!"%identifier)

    return pred, nbres


def writeHeaders(outputprefix, evaluationfile):
    """ Write output file headers """
    outputfilefiltered = "%s.filteredoutput.txt"%outputprefix
    if evaluationfile:
        # read from evaluationfile
        outputfileeval = "%s.filteredoutput.eval.txt"%outputprefix
    else:
        outputfileeval = None

    with open(outputfilefiltered,"a") as o:
        o.write("#identifier diversity     direction viterbiscore indexpred        state  res1  res2\n")
    if evaluationfile:
        with open(outputfileeval,"a") as o:
            o.write("#identifier diversity     direction viterbiscore indexpred        state  res1  res2  T/F   total#rescontacts\n")

    return outputfilefiltered, outputfileeval


def writeOutputBeforeExiting(outputprefix, identifier, diversityvalue, evaluationfile):
    """ In case we run into an error, record something before exiting """
    outputfilefiltered = "%s.filteredoutput.txt"%outputprefix

    # If no Viterbi path has been retained, still record a result line
    with open(outputfilefiltered,"a") as o:
        o.write("%11s %9.2f %13s %12s %9s %12s %5s %5s\n"%(identifier,diversityvalue,"NA","NA","NA","NA","NA","NA"))

    totalnumberrescontacts = "NA"
    if evaluationfile:
        # read from evaluationfile
        outputfileeval = "%s.filteredoutput.eval.txt"%outputprefix
        foundPDB = False
        nextSeq = 0
        evaldic = {} # evaldic = {1: sequence, 2: DSSP, 3: partners1, 4: partners2}
        with open(evaluationfile, "r") as f:
            for l in f:
                if foundPDB:
                    break
                ls = l.split()
                if len(ls)==0 or l[0]=="#":
                    nextSeq = 0
                    continue
                if len(ls) >= 2 and nextSeq > 0:
                    seq = []
                    for t in ls:
                        seq.append(t.strip())
                    evaldic[nextSeq] = seq
                    nextSeq += 1
                    if nextSeq==5:
                        foundPDB = True
                elif len(ls)==1 and (ls[0]==identifier or ls[0]==identifier.upper() or ls[0].replace("_","")==identifier or ls[0].replace("_","")==identifier.upper()):
                    nextSeq = 1
        # register the true contacting pairs
        realresiduecontacts = collections.defaultdict(lambda: collections.defaultdict(lambda: "F"))
        if len(evaldic)==0:
            print("\nWARNING: PDB identifier %s was not found in evaluation file %s, no evaluation will be performed\n"%(identifier, evaluationfile))
            evaluationfile = False
        else:
            for part in [evaldic[3], evaldic[4]]:
                for (ite, p) in enumerate(part):
                    if p=="0":
                        continue
                    realresiduecontacts[ite+1][int(p)] = "T"
                    realresiduecontacts[int(p)][ite+1] = "T"
        # register the total number of contacting pairs
        totalnumberrescontacts = 0
        for i in realresiduecontacts:
            for j in realresiduecontacts[i]:
                if j>=i:
                    continue
                totalnumberrescontacts += 1

        with open(outputfileeval,"a") as o:
            o.write("%11s %9.2f %13s %12s %9s %12s %5s %5s %3s %10i\n"%(identifier,diversityvalue,"NA","NA","NA","NA","NA","NA","NA",totalnumberrescontacts))


def writeOutput(outputfile, identifier, diversityvalue, mostprobablepathendprob, mostprobablepath, indexpredicted):
    """ Write output from path to outputfile """
    with open(outputfile,"a") as o:
        for it,jj in reversed(list(enumerate(mostprobablepath[:-1]))):
            if jj[1]!="start":
                direc = jj[0]
                res1 = int(jj[2])
                res2 = int(jj[3])
                state = jj[1]
                if it == 0 and state=="internal":
                    state = "last"
                elif it == (len(mostprobablepath)-2) and state=="internal":
                    state = "first"
                o.write("%11s %9.2f %13s %12.6f %9i %12s %5i %5i\n"%(identifier,diversityvalue,direc,mostprobablepathendprob,indexpredicted,state,res1,res2))


def writeOutputEval(outputfile, identifier, diversityvalue, mostprobablepathendprob, mostprobablepath, indexpredicted, realresiduecontacts, totalnumberrescontacts):
    """ Write output from path to outputfile with evaluation """
    with open(outputfile,"a") as o:
        for it,jj in reversed(list(enumerate(mostprobablepath[:-1]))):
            if jj[1]!="start":
                direc = jj[0]
                res1 = int(jj[2])
                res2 = int(jj[3])
                state = jj[1]
                if it == 0 and state=="internal":
                    state = "last"
                elif it == (len(mostprobablepath)-2) and state=="internal":
                    state = "first"
                o.write("%11s %9.2f %13s %12.6f %9i %12s %5i %5i %3s %10i\n"%(identifier,diversityvalue,direc,mostprobablepathendprob,indexpredicted,state,res1,res2,realresiduecontacts[res1][res2],totalnumberrescontacts))
