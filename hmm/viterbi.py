""" Viterbi algorithm for HMM decoding """

import collections
import hmm.paramshmm

def initViterbiParallel(viterbi, nbres, dsspMasking, secstructdic, secstructseq, prioroffset, secprior):
    """ Initialize the Viterbi algorithm for the parallel case """
    for i in range(1,nbres+1):
        if (not dsspMasking or i in secstructdic["allprev"]):
            viterbi["start"][i][0] = 0.0 + prioroffset[i] + 0.0 # no secprior for non-assigned secondary structure
    for j in range(1,nbres+1):
        for i in range(j+1,nbres+1):
            # local viterbi
            if not dsspMasking or (i in secstructdic["allprev"] and j in secstructdic["allprev"] and secstructdic[i+1]!=secstructdic[j+1]):
                # For the DSSP case, allow only the combinations where both residues belong to DSSP strands and are not in the same strand
                # For the PSIPRED case, give a start score to every cell
                # The cell (i+1, j+1) will inherit from this start state so the prior offset is calculated at (i+1)-(j+1)=i-j
                # The secondary structure "prior" is calculated at (i,j) but the secstructseq indices start at 0
                # whereas the matrix indices start at 1, hence the secstructseq[i-1]+secstructseq[j-1]
                viterbi["start"][i][j] = 0.0 + prioroffset[i-j] + secprior[secstructseq[i-1]+secstructseq[j-1]]
    return viterbi

def initViterbiAntiparallel(viterbi, nbres, dsspMasking, secstructdic, secstructseq, prioroffset, secprior):
    """ Initialize the Viterbi algorithm for the antiparallel case """
    for i in range(1,nbres+1):
        if (not dsspMasking or (i in secstructdic["allprev"] and i in secstructdic["allprevanti"])):
                viterbi["start"][i][i] = 0.0 + prioroffset[2] + 0.0 # no secprior for non-assigned secondary structure
    for j in range(2,nbres+1): # for j in range(1,nbres+1): # cannot have a start at j=1
        for i in range(j+1,nbres+1):
            # local viterbi
            if not dsspMasking or (i in secstructdic["allprev"] and j in secstructdic["allprevanti"] and secstructdic[i+1]!=secstructdic[j-1]):
                # For the DSSP case, allow only the combinations where both residues belong to DSSP strands and are not in the same strand
                # For the PSIPRED case, give a start score to every cell
                # The cell (i+1, j-1) will inherit from this start state so the prior offset is calculated at (i+1)-(j-1)=i-j+2
                # The secondary structure "prior" is calculated at (i,j) but the secstructseq indices start at 0
                # whereas the matrix indices start at 1, hence the secstructseq[i-1]+secstructseq[j-1]
                viterbi["start"][i][j] = 0.0 + prioroffset[i-j+2] + secprior[secstructseq[i-1]+secstructseq[j-1]]

    return viterbi


def getInitialCrossoutParallel(dsspMasking,secstructdic,nbres,maskarounddiagparallel):
    """ Define the initial set of crossed-out cells for the parallel case """
    crossout = set()
    for i in range(0,nbres+1):
        for jj in range(0,maskarounddiagparallel+1):
            # protecting the region close to the diagonal from getting contacts propagating from a larger distance through bulges
            crossout.add((i+jj,i))
            crossout.add((i,i+jj))
            crossout.add((i-jj,i))
            crossout.add((i,i-jj))
    if dsspMasking:
        # this is just to speed up the Viterbi algorithm,
        # as these regions would get Viterbi score = -inf anyway
        # due to the secondary structure part of the emission probabilities
        for i in range(1,nbres+1):
            for j in range(1,nbres+1):
                if i not in secstructdic["all"] or j not in secstructdic["all"] or (i in secstructdic["all"] and j in secstructdic["all"] and secstructdic[i]==secstructdic[j]):
                    crossout.add((i,j))
    return crossout

def getInitialCrossoutAntiparallel(dsspMasking,secstructdic,nbres,maskarounddiagantiparallel):
    """ Define the initial set of crossed-out cells for the antiparallel case """
    crossout = set()
    for i in range(0,nbres+1):
        for jj in range(0,maskarounddiagantiparallel+1):
            # protecting the region close to the diagonal from getting contacts propagating from a larger distance through bulges
            crossout.add((i+jj,i))
            crossout.add((i,i+jj))
            crossout.add((i-jj,i))
            crossout.add((i,i-jj))

    if dsspMasking:
        # this is just to speed up the Viterbi algorithm,
        # as these regions would get Viterbi score = -inf anyway
        # due to the secondary structure part of the emission probabilities
        for i in range(1,nbres+1):
            for j in range(1,nbres+1):
                if i not in secstructdic["all"] or j not in secstructdic["all"] or (i in secstructdic["all"] and j in secstructdic["all"] and secstructdic[i]==secstructdic[j]):
                    crossout.add((i,j))
    return crossout


def getCrossoutFromAllpaths(direction, crossout, newcrossout, ite, allpaths, nbres, numcontacts, maskaroundcontact):
    """ Get the current set of crossed-out cells
    depending on the predictions that we have already made (contained in allpaths) """
    # Crossout a region of width maskaroundcontact on each side of the main EE contact
    mask = maskaroundcontact
    for (i,(prob,path)) in enumerate(allpaths):
        if i+2<ite:
            # we have seen this path before, so the cells have been crossed-out already
            continue
        for k in path[:-1]: # the last path state is the start state
            direc = k[0]
            if direc!=direction:
                # we only crossout the area around this contact for the corresponding directions, so that
                # we have both the parallel and antiparallel contacts in the final (non-filtered) list of best Viterbi paths
                # if both have a Viterbi score above the threshold
                continue
            state = k[1]
            crossout.add((k[2],k[3]))
            # newcrossout contains only the cells that have been crossed-out while running the Viterbi algorithm,
            # not the ones that were initially crossed-out (e.g. close to the diagonal or in DSSP non-E regions)
            # (newcrossout is used in prediction-shortening mode, to speed up the recalculation)
            newcrossout.add((k[2],k[3]))
            for ii in range(mask):
                for jj in range(mask):
                    ## if direc=="Antiparallel":
                    crossout.add((k[2]-ii,k[3]-jj))
                    crossout.add((k[2]+ii,k[3]+jj))
                    newcrossout.add((k[2]-ii,k[3]-jj))
                    newcrossout.add((k[2]+ii,k[3]+jj))
                    ## elif direc=="Parallel":
                    crossout.add((k[2]-ii,k[3]+jj))
                    crossout.add((k[2]+ii,k[3]-jj))
                    newcrossout.add((k[2]-ii,k[3]+jj))
                    newcrossout.add((k[2]+ii,k[3]-jj))
    return crossout, numcontacts, newcrossout


def updateViterbiCrossout(dsspMasking, direction, crossout, viterbi, viterbiptr, i, j, trprob, secstructseq, nbres, emissions, recalculateCount):
    """ Update the Viterbi scores and pointers for cell (i,j) """
    if direction=="Parallel":
        jprev = j-1
    elif direction=="Antiparallel":
        jprev = j+1

    # for states where we take a step in both directions
    for state in ["bulgei0","bulgej0","first","internal","last"]:
        # stprev = "start" is the extra element over which we maximize that enables the contact (HMM path) to start here at (i,j) (local Viterbi)
        mylist = [(viterbi[stprev][i-1][jprev] + (trprob[stprev][state]),stprev) for stprev in ["bulgei0","bulgei1","bulgei2","bulgej0","bulgej1","bulgej2","start","first","internal","last"] if trprob[stprev][state] != (-float("inf"))]
        (mymaxscore,mymaxst) = max(mylist)
        viterbi[state][i][j] = mymaxscore + emissions[direction][i][j][state]
        viterbiptr[state][i][j] = (mymaxst,i-1,jprev)

    # for states where we take a step in only one direction
    # bulgei1 can only derive from bulgei0 (and trprob is 1)
    state="bulgei1"
    viterbi[state][i][j] = viterbi["bulgei0"][i][jprev] + (emissions[direction][i][j][state])
    viterbiptr[state][i][j] = ("bulgei0",i,jprev)
    # bulgej1 can only derive from bulgej0 (and trprob is 1)
    state="bulgej1"
    viterbi[state][i][j] = viterbi["bulgej0"][i-1][j] + (emissions[direction][i][j][state])
    viterbiptr[state][i][j] = ("bulgej0",i-1,j)
    # bulgei2 can only derive from bulgei1
    state="bulgei2"
    viterbi[state][i][j] = viterbi["bulgei1"][i][jprev] + (emissions[direction][i][j][state]) + (trprob["bulgei1"][state])
    viterbiptr[state][i][j] = ("bulgei1",i,jprev)
    # bulgej2 can only derive from bulgej1
    state="bulgej2"
    viterbi[state][i][j] = viterbi["bulgej1"][i-1][j] + (emissions[direction][i][j][state]) + (trprob["bulgej1"][state])
    viterbiptr[state][i][j] = ("bulgej1",i-1,j)

    return viterbi,viterbiptr


def getNextViterbi(viterbi, viterbiptr, viterbirecalclower, crossout, myitems, itestart, problematic):
    """ Get the next best Viterbi score which is compatible with crossed-out regions and retrieve the associated path """
    if itestart >= len(myitems):
        # we have reached the end of myitems
        return [], -float("inf"), itestart, False

    ite = itestart
    containsCrossout = True
    probEncountered = False

    while len(myitems)>=1 and ite < len(myitems) and containsCrossout:
        # Loop over myitems, starting at position itestart, until we reach a path that does not contain crossed-out pairs
        containsCrossout = False
        mymax2 = myitems[ite]

        ite += 1
        (p,dirfinal,stfinal,ifinal,jfinal) = mymax2
        if (ifinal,jfinal) in problematic[dirfinal]:
            # (used for prediction-shortening mode): we have encountered a problematic pair of residues (i.e. a pair that belonged to a long prediction)
            # so the Viterbi scores should be recalculated (otherwise we risk to discard interesting paths)
            probEncountered = True

        if (ifinal,jfinal) in crossout[dirfinal]:
            # this pair was crossed-out before, discard this path
            containsCrossout = True
            continue

        # Back-trace, looking for crossed-out pairs which lead to exclusion of the path
        dircurr = dirfinal
        stcurr = stfinal
        icurr = ifinal
        jcurr = jfinal
        while stcurr != "start":
            (stcurr,icurr,jcurr) = viterbiptr[dircurr][stcurr][icurr][jcurr]

            if (icurr,jcurr) in problematic[dirfinal]:
                # (used for prediction-shortening mode): we have encountered a problematic pair of residues (i.e. a pair that belonged to a long prediction)
                # so the Viterbi scores should be recalculated (otherwise we risk to discard interesting paths)
                probEncountered = True

            if stcurr != "start" and (icurr,jcurr) in crossout[dirfinal]:
                (stcurr,icurr,jcurr) = viterbiptr[dircurr][stcurr][icurr][jcurr]
                # this pair was crossed-out before, discard this path
                containsCrossout = True
                # break here, this saves some time on the backtracing
                break

    if not containsCrossout:
        mymax = mymax2
    else:
        # we have reached the end of myitems without finding a path not containing crossed-out pairs
        return [], -float("inf"), itestart, False

    if mymax[0] <= viterbirecalclower:
        # we have reached the minimum Viterbi score that we set as a threshold
        return [], -float("inf"), itestart, False

    # Back-trace and register the path
    (p,dirfinal,stfinal,ifinal,jfinal) = mymax
    mostprobablepath = [[dirfinal,stfinal,ifinal,jfinal]]
    mostprobablepathendprob = p
    dircurr = dirfinal
    stcurr = stfinal
    icurr = ifinal
    jcurr = jfinal
    while stcurr != "start":
        (stcurr,icurr,jcurr) = viterbiptr[dircurr][stcurr][icurr][jcurr]
        mostprobablepath.append([dircurr,stcurr,icurr,jcurr])

    return  mostprobablepath,mostprobablepathendprob, ite, probEncountered


def runViterbi(nbres, trprob, prioroffset, secprior, secprobdic, emissions, viterbirecalclower, viterbirecalclowerPSM, PSMparams, dsspMasking, secstructseq, secstructdic, maskarounddiagparallel, maskarounddiagantiparallel, maskaroundcontact, noshorteningmode):
    recalculateViterbi = True
    recalculateCount = 0
    recalculateFactor = 0
    problematic = collections.defaultdict(set) # residue pairs contained in the paths that are too long
    fullliststates = ["bulgei0","bulgei1","bulgei2","bulgej0","bulgej1","bulgej2","internal","start","end","first","last"]

    # Loop until we do not have to run the Viterbi algorithm anymore
    # (the loop will be traversed only once if we do not enter the prediction-shortening mode)
    while recalculateViterbi:
        recalculateViterbi = False # will only become True again if we enter the prediction-shortening mode
        ite = 0
        maxite = (len(fullliststates)*nbres*(nbres-1))/2 # number of values in the viterbi dictionary
        p = viterbirecalclower + 1

        # if we are recalculating, update the transition probabilities
        if recalculateCount >= 1:
            print("(Prediction-shortening mode) Rerunning the Viterbi algorithm with lower transition probabilities", recalculateCount)
            trprob = hmm.paramshmm.updateTransitionProbabilities(trprob, recalculateFactor, PSMparams["PSMtrprobstepsize"])

        recalculateFactor = 0
        allpaths = []
        crossoutparallel = getInitialCrossoutParallel(dsspMasking, secstructdic, nbres, maskarounddiagparallel)
        crossoutantiparallel = getInitialCrossoutAntiparallel(dsspMasking, secstructdic, nbres, maskarounddiagantiparallel)

        numcontactsparallel = collections.defaultdict(int)
        numcontactsantiparallel = collections.defaultdict(int)

        # Initialization and recursion steps of the Viterbi algorithm (separately for the parallel and the antiparallel case)
        viterbi = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: -float("inf")))))
        viterbiptr = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: ("start",0,0)))))

        direction = "Parallel"
        newcrossoutparallel = set()
        crossoutparallel, numcontactsparallel, newcrossoutparallel = getCrossoutFromAllpaths(direction, crossoutparallel, newcrossoutparallel, ite, allpaths, nbres, numcontactsparallel, maskaroundcontact)
        (viterbi[direction]) = initViterbiParallel(viterbi[direction], nbres, dsspMasking, secstructdic, secstructseq, prioroffset[direction], secprior[direction])
        for j in range(1,nbres+1):
            for i in range(j+1,nbres+1):
                if (i,j) in crossoutparallel:
                    continue
                (viterbi[direction],viterbiptr[direction]) = updateViterbiCrossout(dsspMasking,direction,crossoutparallel,viterbi[direction],viterbiptr[direction],i,j,trprob[direction],secstructseq,nbres,emissions,recalculateCount)

        direction = "Antiparallel"
        newcrossoutantiparallel = set()
        crossoutantiparallel, numcontactsantiparallel, newcrossoutantiparallel = getCrossoutFromAllpaths(direction, crossoutantiparallel, newcrossoutantiparallel, ite, allpaths, nbres, numcontactsantiparallel, maskaroundcontact)
        (viterbi[direction]) = initViterbiAntiparallel(viterbi[direction], nbres, dsspMasking, secstructdic, secstructseq, prioroffset[direction], secprior[direction])
        for k in range(1,nbres):
            # k is the difference between i and j
            for j in range(1,nbres-k+1):
                i = j+k
                if (i,j) in crossoutantiparallel:
                    continue
                (viterbi[direction],viterbiptr[direction]) = updateViterbiCrossout(dsspMasking,direction,crossoutantiparallel,viterbi[direction],viterbiptr[direction],i,j,trprob[direction],secstructseq,nbres,emissions,recalculateCount)

        # Termination of the Viterbi algorithm
        # To speed-up the termination, build a list of states that have non-zero transition probabilities towards the "end" state
        liststatescanend = {}
        for k in ["Parallel","Antiparallel"]:
            liststatescanend[k] = []
            for k1 in fullliststates:
                if trprob[k][k1]["end"] > (-float("inf")):
                    liststatescanend[k].append(k1)

        # Crossout for both parallel and antiparallel cases.
        crossout = {}
        crossout["Parallel"] = crossoutparallel
        crossout["Antiparallel"] = crossoutantiparallel

        multi = {"Parallel":1, "Antiparallel":(-1)}

        # This list merges the parallel and antiparallel cases and calculates on-the-fly the final Viterbi scores including the termination step.
        myitems = [(viterbi[k][k1][k2][k3]+(trprob[k][k1]["end"])+  (secprobdic[k]["end"][secstructseq[k2-1+1]+secstructseq[k3-1+multi[k]*1]+secstructseq[k2-1]+secstructseq[k3-1]]), k, k1, k2, k3) for k in viterbi for k1 in liststatescanend[k] for k2 in viterbi[k][k1] for k3 in viterbi[k][k1][k2] if ((viterbi[k][k1][k2][k3]+1 > viterbirecalclower))]

        # Sort in decreasing order, so the best Viterbi score will correspond to the first element of myitems.
        myitems.sort(reverse=True)

        itestart = 0

        # Loop over elements of myitems and retrieve the successive best Viterbi paths which are compatible with previous predictions
        while ite < maxite and p > viterbirecalclower:
            ite +=1

            # Update the crossout dictionary
            crossoutparallel, numcontactsparallel, newcrossoutparallel = getCrossoutFromAllpaths("Parallel",crossoutparallel, newcrossoutparallel,ite,allpaths,nbres,numcontactsparallel, maskaroundcontact)
            crossoutantiparallel, numcontactsantiparallel, newcrossoutantiparallel = getCrossoutFromAllpaths("Antiparallel",crossoutantiparallel, newcrossoutantiparallel,ite,allpaths,nbres,numcontactsantiparallel, maskaroundcontact)

            crossout = {}
            crossout["Parallel"] = crossoutparallel
            crossout["Antiparallel"] = crossoutantiparallel

            if len(myitems)==0:
                # reached end of myitems list
                break

            # Get the next best Viterbi path
            (path,pathprob,itestart,probEncountered) = getNextViterbi(viterbi, viterbiptr, viterbirecalclower-PSMparams["PSMtrprobstepsize"]*recalculateCount*10, crossout, myitems, itestart, problematic)

            ### Start of PREDICTION-SHORTENING MODE (PSM)-specific code (rerun the Viterbi algorithm) ###
            # Only if recalculateCount is >= 1 AND probEncountered is True do we rerun the Viterbi algorithm
            if ite > 1 and recalculateCount >= 1 and probEncountered: # if ite==1, we have calculated the Viterbi scores just before entering the while loop
                # set viterbirecalclower to something higher to save time
                if viterbirecalclowerPSM > viterbirecalclower:
                    viterbirecalclower = viterbirecalclowerPSM
                # also set the maximum number of meaningful predictions to 100 to save time
                maxite = 100

                # Store the previous Viterbi dictionary, some of the Viterbi scores will not be recalculated
                oldviterbi = viterbi
                oldviterbiptr = viterbiptr
                # Initialize the Viterbi dictionaries
                viterbi = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: -float("inf")))))
                viterbiptr = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: ("start",0,0)))))

                # Initialization and recursion of the Viterbi algorithm for the parallel case
                direction = "Parallel"
                (viterbi[direction]) = initViterbiParallel(viterbi[direction], nbres, dsspMasking, secstructdic, secstructseq, prioroffset[direction], secprior[direction])
                crossoutSeen = False
                for j in range(1,nbres+1):
                    for i in range(j+1,nbres+1):
                        if (i,j) in newcrossoutparallel:
                            if not crossoutSeen:
                                # find the first crossed-out position to determine which Viterbi scores should be recalculated
                                (i0,j0) = (i,j)
                            crossoutSeen = True
                            continue
                        if (i,j) in crossoutparallel:
                            continue
                        if not crossoutSeen or i < i0:
                            # for cells not impacted by the cross-out, do not recalculate the Viterbi scores
                            for state in fullliststates:
                                viterbi[direction][state][i][j] = oldviterbi[direction][state][i][j]
                                viterbiptr[direction][state][i][j] = oldviterbiptr[direction][state][i][j]
                        else:
                            (viterbi[direction],viterbiptr[direction]) = updateViterbiCrossout(dsspMasking,direction,crossoutparallel,viterbi[direction],viterbiptr[direction],i,j,trprob[direction],secstructseq,nbres,emissions,recalculateCount)

                # Initialization and recursion of the Viterbi algorithm for the antiparallel case
                direction = "Antiparallel"
                (viterbi[direction]) = initViterbiAntiparallel(viterbi[direction], nbres, dsspMasking, secstructdic, secstructseq, prioroffset[direction], secprior[direction])
                crossoutSeen = False
                for k in range(1,nbres):
                    # k is the difference between i and j
                    for j in range(1,nbres-k+1):
                        i = j+k
                        if (i,j) in newcrossoutantiparallel:
                            if not crossoutSeen:
                                # find the first crossed-out position to determine which Viterbi scores should be recalculated
                                (i0,j0) = (i,j)
                            crossoutSeen = True
                            continue
                        if (i,j) in newcrossoutantiparallel:
                            continue
                        if (not crossoutSeen) or (i < i0 or j > j0):
                            # for cells not impacted by the cross-out, do not recalculate the Viterbi scores
                            for state in fullliststates:
                                viterbi[direction][state][i][j] = oldviterbi[direction][state][i][j]
                                viterbiptr[direction][state][i][j] = oldviterbiptr[direction][state][i][j]
                        else:
                            (viterbi[direction],viterbiptr[direction]) = updateViterbiCrossout(dsspMasking,direction,crossoutantiparallel,viterbi[direction],viterbiptr[direction],i,j,trprob[direction],secstructseq,nbres,emissions,recalculateCount)

                # Termination of the Viterbi algorithm
                # To speed-up the termination, build a list of states that have non-zero transition probabilities towards the "end" state
                liststatescanend = {}
                for k in ["Parallel","Antiparallel"]:
                    liststatescanend[k] = []
                    for k1 in fullliststates:
                        if trprob[k][k1]["end"] > (-float("inf")):
                            liststatescanend[k].append(k1)

                # Crossout for both parallel and antiparallel cases.
                crossout = {}
                crossout["Parallel"] = crossoutparallel
                crossout["Antiparallel"] = crossoutantiparallel

                multi = {"Parallel":1, "Antiparallel":(-1)}

                # This list merges the parallel and antiparallel cases and calculates on-the-fly the final Viterbi scores including the termination step.
                myitems = [(viterbi[k][k1][k2][k3]+(trprob[k][k1]["end"])+ (secprobdic[k]["end"][secstructseq[k2-1+1]+secstructseq[k3-1+multi[k]*1]+secstructseq[k2-1]+secstructseq[k3-1]]), k, k1, k2, k3) for k in viterbi for k1 in liststatescanend[k] for k2 in viterbi[k][k1] for k3 in viterbi[k][k1][k2] if ((viterbi[k][k1][k2][k3]+1+PSMparams["PSMtrprobstepsize"]*recalculateCount*10 > viterbirecalclower))]
                # Sort in decreasing order, so the best Viterbi score will correspond to the first element of myitems.
                myitems.sort(reverse=True)

                itestart = 0 # reset the itestart because we have recalculated some of the Viterbi scores
                if len(myitems)==0:
                    break

                # Get the next best Viterbi path
                (path,pathprob,itestart, probEncountered) = getNextViterbi(viterbi, viterbiptr, viterbirecalclower-PSMparams["PSMtrprobstepsize"]*recalculateCount*10, crossout, myitems, itestart, problematic)
            ### End of PREDICTION-SHORTENING MODE (PSM)-specific code (rerun the Viterbi algorithm) ###

            # Should we enter prediction-shortening mode?
            if not noshorteningmode and not dsspMasking and ite < PSMparams["PSMmaxitenumber"] and len(path) > 1 and ((path[::-1][1][0]=="Parallel" and len(path) > PSMparams["PSMlongestparallel"]) or (path[::-1][1][0]=="Antiparallel" and len(path) > PSMparams["PSMlongestantiparallel"])):
                for k in path[::-1][1:]:
                    res1 = k[2]
                    res2 = k[3]
                    problematic[k[0]].add((res1,res2))
                recalculateViterbi = True
                if recalculateCount == 0:
                    print("\nWARNING: Now entering prediction-shortening mode, the prediction might take a few minutes to several hours (depending on the protein size).\nYou might consider switching this mode off with option -l (at the potential cost of more false positives due to very long predicted contacts.\n")
                if len(path) > PSMparams["PSMspeeduppathlength"]:
                    #speed up process by substracting twice the PSMtrprobstepsize from the transition probs
                    recalculateFactor = 1
                break

            # Re-shift the Viterbi score for the path so that the results are comparable with results from the "regular" version of the script.
            # With the "regular" version, recalculateCount is 0 so this will not change the Viterbi score.
            pathprob += (PSMparams["PSMtrprobstepsize"]) * recalculateCount * (len(path)-2) # Add PSMtrprobstepsize for each transition except the transition to the last state
            allpaths.append((pathprob,path))
            if len(path) > 1:
                print(ite,"%13.10f"%pathprob,path[::-1][1][0],)
                for k in path[::-1][1:]:
                    res1 = k[2]
                    res2 = k[3]
                    print(k[1],k[2],k[3],"-",)
                print("")

            p = pathprob

        # increase the recalculateCount counter
        recalculateCount += (1 + recalculateFactor)

    if recalculateCount > 1:
        # we need to sort the allpaths table because we fiddled with the scores
        allpaths.sort(reverse=True)

    return allpaths, recalculateCount

