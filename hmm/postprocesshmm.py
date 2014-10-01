""" Post-processing HMM decoding results and writing output """

import collections
import hmm.iohmm

def postProcessPaths(allpaths, outputprefix, minimumprob, identifier, diversityvalue, mask, secstructdic, evaluationfile):
    """ Post process the list of all retained Viterbi paths and write output """
    # Write headers
    outputfilefiltered, outputfileeval = hmm.iohmm.writeHeaders(outputprefix, evaluationfile)

    # If no Viterbi path has been retained, still record a result line
    if len(allpaths)==0 or allpaths == [(-float("inf"), [])]:
        hmm.iohmm.writeOutputBeforeExiting(outputprefix, identifier, diversityvalue, evaluationfile)
        
    # Read real contacts from evaluation file
    totalnumberrescontacts = "NA"
    if evaluationfile:
        # read from evaluationfile
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
            print "\nWARNING: PDB identifier %s was not found in evaluation file %s, no evaluation will be performed\n"%(identifier, evaluationfile)
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

    # Run the post-processing
    predictedcontacts = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(float)))
    predictedcontactpaths = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))
    indexpredicted = 0
    numcontactsperres = collections.defaultdict(int)
    predictedpairs = set()
    filteringcrossout = set()
    predictedstrandpairs = set()

    for (mostprobablepathendprob,mostprobablepath) in allpaths:
        if mostprobablepathendprob <= minimumprob:
            continue

        if len(mostprobablepath)==1:
            # the path only contains a "start" state
            continue

        ressetforcontacts = set()
        currentpairs = set()
        currentstrandpairs = set()

        direc = mostprobablepath[-1][0]
        indexpredicted += 1

        # Loop through residues involved in this path
        for it,jj in reversed(list(enumerate(mostprobablepath[:-1]))):
            if jj[1]!="start":
                res1 = int(jj[2])
                res2 = int(jj[3])
                # Record some data for the filtering step
                ressetforcontacts.add(res1)
                ressetforcontacts.add(res2)
                currentpairs.add((res1,res2))
                if res1 in secstructdic and res2 in secstructdic:
                    currentstrandpairs.add((secstructdic[res1], secstructdic[res2]))
                    currentstrandpairs.add((secstructdic[res2], secstructdic[res1]))

        # Filtering: for DSSP, cannot predict an interaction between two residues belonging to the same strand
        if len(currentstrandpairs) !=0 and secstructdic[res1]==secstructdic[res2]:
            continue
        # Filtering: for DSSP, cannot predict the interaction between 2 strands twice
        shouldContinue = False
        for p in currentstrandpairs:
            if p in predictedstrandpairs:
                shouldContinue = True
                break
        if shouldContinue:
            continue

        # Filtering: maximum 2 contacts per residue
        shouldContinue = False
        for r in ressetforcontacts:
            if numcontactsperres[r]>=2:
                shouldContinue = True
                break
        if shouldContinue:
            continue

        # Filtering: cannot predict the same pair in two different strand-strand contacts
        for (r1,r2) in currentpairs:
            if (r1,r2) in predictedpairs or (r2,r1) in predictedpairs:
                shouldContinue = True
                break
            if (r1,r2) in filteringcrossout or (r2,r1) in filteringcrossout:
                shouldContinue = True
                break
        if shouldContinue:
            continue

        # Filtering: once this contact has been accepted (given what was seen before),
        # update the number of contacts and the predictedpairs and the crossout
        for r in ressetforcontacts:
            numcontactsperres[r] += 1
        for p in currentpairs:
            predictedpairs.add(p)
        for p in currentstrandpairs:
            predictedstrandpairs.add(p)
        (istart,jstart) = (mostprobablepath[-1][2],mostprobablepath[-1][3])
        direc = mostprobablepath[0][0]
        if direc != "Antiparallel":
            (iend,jend) = (mostprobablepath[0][2]+1,mostprobablepath[0][3]+1)
        else:
            (iend,jend) = (mostprobablepath[0][2]+1,mostprobablepath[0][3]-1)
        for ki in range(istart,iend+1):
            for kj in range(min(jstart,jend),max(jstart,jend)+1):
                for ii in xrange(mask):
                    for jj in xrange(mask):
                        if direc=="Antiparallel":
                            if (ki-ii) in range(istart,iend+1) and (kj-jj) in range(min(jstart,jend),max(jstart,jend)+1):
                                filteringcrossout.add((ki-ii,kj-jj))
                            if (ki+ii) in range(istart,iend+1) and (kj+jj) in range(min(jstart,jend),max(jstart,jend)+1):
                                filteringcrossout.add((ki+ii,kj+jj))
                        elif direc=="Parallel":
                            if (ki-ii) in range(istart,iend+1) and (kj+jj) in range(min(jstart,jend),max(jstart,jend)+1):
                                filteringcrossout.add((ki-ii,kj+jj))
                            if (ki+ii) in range(istart,iend+1) and (kj-jj) in range(min(jstart,jend),max(jstart,jend)+1):
                                filteringcrossout.add((ki+ii,kj-jj))


        # Filtering: write filtered output (only accepted contacts)
        hmm.iohmm.writeOutput(outputfilefiltered, identifier, diversityvalue, mostprobablepathendprob, mostprobablepath, indexpredicted)
        if evaluationfile:
            # Write filtered output with evaluation results (residue-level)
            hmm.iohmm.writeOutputEval(outputfileeval, identifier, diversityvalue, mostprobablepathendprob, mostprobablepath, indexpredicted, realresiduecontacts, totalnumberrescontacts)
