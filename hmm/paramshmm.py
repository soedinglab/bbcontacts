""" Retrieve and adjust parameters for the HMM decoding """

import collections
import numpy as np

# Do not print warnings when calculating np.log(0)
np.seterr(divide='ignore')

from scipy.special import gamma


def smoothMatrix(predcouplings, nbres, smoothingsize):
    """ Smooth the coupling matrix with the given SMOOTHINGSIZE parameter """
    localmean =  collections.defaultdict(lambda: collections.defaultdict(float))
    for i in range(1, nbres+1):
        for j in range(1, nbres+1):
            if j<=i:
                continue
            localsum = [] # intermediate list to store values retrieved from the smoothing area
            for kk in range(0, smoothingsize+1):
                for ll in range(0, smoothingsize+1):
                    try:
                        localsum.append(predcouplings[i+kk][j+ll])
                    except KeyError:
                        pass
                    if ll != 0: # otherwise some couplings count twice
                        try:
                            localsum.append(predcouplings[i+kk][j-ll])
                        except KeyError:
                            pass
                    if kk == 0:
                        continue # otherwise some couplings count twice
                    try:
                        localsum.append(predcouplings[i-kk][j+ll])
                    except KeyError:
                        pass
                    if ll != 0: # otherwise some couplings count twice
                        try:
                            localsum.append(predcouplings[i-kk][j-ll])
                        except KeyError:
                            pass
            # localmean contains the average over the area centered at (i,j) and extending by SMOOTHINGSIZE in each direction
            localmean[i][j] = sum(localsum) / float(len(localsum))
            localmean[j][i] = localmean[i][j]
    for i in range(1, nbres+1):
        for j in range(1, nbres+1):
            # subtract the localmean from each cell of the coupling matrix
            predcouplings[i][j] -= localmean[i][j]
    return predcouplings


def getSecondaryStructureProbabilities(inputfile):
    """ Retrieve secondary structure probabilities """
    sec = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(float)))
    half = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(float))))
    with open(inputfile,"r") as f:
        for l in f:
            if l[0]=="#":
                continue
            line = l.split()
            direction = line[0]
            state = line[3]
            t = line[1] # current sec struct state
            tp = line[2] # (given) previous sec struct state
            half[direction][state][t][tp] = float(l.split()[4])

    for d in half.keys():
        for s in half[d].keys():
            if len(s)>=6 and s[0:6] in ["bulgei1","bulgei2"]:
                for ti in ["C","E","H"]:
                    for tj in ["C","E","H"]:
                        for tip in ["C","E","H"]:
                            for tjp in ["C","E","H"]:
                                # for bulge extension states (buylgex1, bulgex2) count only the part of the sec struct corresponding to the long side of the bulge
                                sec[d][s][ti+tj+tip+tjp] = half[d][s][tj][tjp]
            elif len(s)>=6 and s[0:6] in ["bulgej1","bulgej2"]:
                for ti in ["C","E","H"]:
                    for tj in ["C","E","H"]:
                        for tip in ["C","E","H"]:
                            for tjp in ["C","E","H"]:
                                # for bulge extension states (buylgex1, bulgex2) count only the part of the sec struct corresponding to the long side of the bulge
                                sec[d][s][ti+tj+tip+tjp] = half[d][s][ti][tip]
            else:
                for ti in ["C","E","H"]:
                    for tj in ["C","E","H"]:
                        for tip in ["C","E","H"]:
                            for tjp in ["C","E","H"]:
                                # for other states, take the product of the conditional probs
                                sec[d][s][ti+tj+tip+tjp] = half[d][s][ti][tip] * half[d][s][tj][tjp]
            if s=="background":
                for ti in ["C","E","H"]:
                    for tj in ["C","E","H"]:
                        for tip in ["C","E","H"]:
                            for tjp in ["C","E","H"]:
                                # define a special background for bulgex1 and bulgex2 states which also only counts one part of the probability
                                sec[d]["backgroundforbulgei"][ti+tj+tip+tjp] = half[d][s][tj][tjp]
                                sec[d]["backgroundforbulgej"][ti+tj+tip+tjp] = half[d][s][ti][tip]
        # Calculate once and for all the log-odds with respect to the background
        for s in half[d].keys():
            if s in ["backgroundforbulgei","backgroundforbulgej","background"]:
                continue
            for ti in ["C","E","H"]:
                for tj in ["C","E","H"]:
                    for tip in ["C","E","H"]:
                        for tjp in ["C","E","H"]:
                            if s in ["bulgei1","bulgei2"]:
                                sec[d][s][ti+tj+tip+tjp] = np.log(sec[d][s][ti+tj+tip+tjp] / sec[d]["backgroundforbulgei"][ti+tj+tip+tjp])
                            elif s in ["bulgej1","bulgej2"]:
                                sec[d][s][ti+tj+tip+tjp] = np.log(sec[d][s][ti+tj+tip+tjp] / sec[d]["backgroundforbulgej"][ti+tj+tip+tjp])
                            else:
                                sec[d][s][ti+tj+tip+tjp] = np.log(sec[d][s][ti+tj+tip+tjp] / sec[d]["background"][ti+tj+tip+tjp])

    return sec

def getSecondaryStructurePrior(inputfile):
    """ Retrieve non-conditional secondary structure probabilities that will be used for the "start" state """
    secprior = collections.defaultdict(lambda: collections.defaultdict(float))
    bgprior = collections.defaultdict(lambda: collections.defaultdict(float))
    with open(inputfile, "r") as f:
        for l in f:
            if l[0]=="#":
                direction = l.split()[1]
                continue
            if l.split()[1]=="background":
                # background always comes first
                bgprior[direction][l.split()[2]] = float(l.split()[3])
            elif l.split()[1]=="start":
                # calculate the log-odds with respect to the background
                secprior[direction][l.split()[2]] = np.log(float(l.split()[3])/bgprior[direction][l.split()[2]])
    return secprior


def getTransformedGammaParameters(pdb, diversityvalue, inputfile):
    """ Retrieve the parameters for the coupling distribution fits
    and make sure they are in the expected range
    (otherwise set lowDiversity to True, so the couplings will not count in the HMM) """

    params = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(float))))

    # We set a minimum value for the parameters which should not get negative
    minimumflist = ["alpha_pos","beta_pos","b_pos","alpha_neg","b_neg","w"]
    minimumf = 0.0
    usefulm = [1,2,3]

    lowDiversity = False

    liststates = ["bulgei0","bulgei1","bulgei2","bulgej0","bulgej1","bulgej2","first","internal","last"]
    paramlist = ["alpha_pos","beta_pos","b_pos","alpha_neg","b_neg","w","shift"]
    directionlist = ["Parallel","Parallel","Parallel","Antiparallel","Antiparallel","Antiparallel"]
    mlist = [["b"],[2],[1,3],["b"],[2],[1,3]]
    with open(inputfile,"r") as f:
        ii = 0
        for l in f:
            if l[0]=="#":
                continue
            if l.split()[0]=="\"alpha_pos_intercept\"":
                # header line
                continue
            if ii >= len(directionlist):
                # stop reading inputfile after len(directionlist) lines
                break
            direction = directionlist[ii]
            for m in mlist[ii]:
                if m=="b":
                    for mm in usefulm: # the background parameters apply to all 3 pattern positions
                        for ite,p in enumerate(paramlist):
                            if p!="shift":
                                # linear fits for all parameters but the shift
                                f = float(l.split()[ite])+ diversityvalue * float(l.split()[ite+6])
                                if p not in minimumflist or f > minimumf:
                                    params[direction]["background"][mm][p] = f
                                else:
                                    lowDiversity = True
                            else:
                                # quadratic fit for the shift
                                f = float(l.split()[12])+ diversityvalue * float(l.split()[13])+ diversityvalue*diversityvalue * float(l.split()[14])
                                params[direction]["background"][mm][p] = f
                else:
                    s = "internal"
                    # the positive parameters are the same for all states,
                    # but different for pattern position 2 (central contact, stronger) and positions 1-3 (side contacts, weaker)
                    for ite,p in enumerate(paramlist):
                        if p!="shift":
                            # linear fits for all parameters but the shift
                            f = float(l.split()[ite])+ diversityvalue * float(l.split()[ite+6])
                            if p not in minimumflist or f > minimumf:
                                params[direction][s][m][p] = f
                            else:
                                lowDiversity = True
                        else:
                            # quadratic fit for the shift
                            f = float(l.split()[12])+ diversityvalue * float(l.split()[13])+ diversityvalue*diversityvalue * float(l.split()[14])
                            params[direction][s][m][p] = f
            ii += 1
    return params, lowDiversity


def getTransitionProbabilities(trprobfile):
    """ Retrieve the transition probabilities """
    d = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: -float("inf"))))
    with open(trprobfile, "r") as f:
        for l in f:
            if l[0] == "#":
                continue
            # d[direction][statestart][stateend] = trprob
            d[l.split()[0]][l.split()[1]][l.split()[2]] = np.log(float(l.split()[3]))
    return d

def updateTransitionProbabilities(trprob, recalculateFactor, PSMtrprobstepsize):
    """ In prediction-shortening mode (PSM):
    each time we encounter a long prediction, the transition probabilities are reduced by PSMtrprobstepsize
    except for the transition probabilities to the "last" state which are just increased
    so that the sum of all transition probabilities starting from any state is always 1"""
    for d in trprob:
        for k in trprob[d]:
            if k not in ["internal","bulgei1","bulgej1","bulgei2","bulgej2","first","last"]:
                continue
            otherprob = 0.0
            for l in trprob[d][k]:
                if l!="last":
                    if l=="end":
                        otherprob += np.exp(trprob[d][k][l])
                    else:
                        trprob[d][k][l] -= PSMtrprobstepsize * (1+recalculateFactor)
                        otherprob += np.exp(trprob[d][k][l])
            trprob[d][k]["last"] = np.log(1.0-otherprob)
    return trprob


def getPriorOffset(prioroffsetfile, nbres):
    """ Retrieve prior offset depending on distance to diagonal """
    prioroffset = collections.defaultdict(lambda: collections.defaultdict(float))
    D1 = collections.defaultdict(int)
    D2 = collections.defaultdict(int)
    linpars = collections.defaultdict(list)
    otherpars = collections.defaultdict(list)
    with open(prioroffsetfile, "r") as f:
        for l in f:
            if l[0] == "#":
                continue
            line = l.split()
            direction = line[0]
            if line[1] == "D1":
                D1[direction] = int(line[2])
                continue
            if line[1] == "D2":
                D2[direction] = int(line[2])
                continue
            if int(line[1]) <= D1[direction]:
                # explicit prior offset up to a distance of D1[direction]
                prioroffset[direction][int(line[1])] = np.log(float(line[2]))
            elif int(line[1]) <= D1[direction]+2:
                # two parameters for the linear dependency region between D1[direction]+1 and D2[direction]
                linpars[direction].append(float(line[2]))
            else:
                # three parameters for the exponential decay after D2[direction]
                otherpars[direction].append(float(line[2]))

    for direction in prioroffset.keys():
        for d in range(D1[direction]+1,D2[direction]+1):
            prioroffset[direction][d] = np.log(linpars[direction][0] + float(d) * linpars[direction][1])
        for d in range(D2[direction]+1,nbres+1):
            prioroffset[direction][d] = np.log(otherpars[direction][0] + otherpars[direction][1] * np.exp(-float(d)/otherpars[direction][2]))

    return prioroffset


def setViterbiThresholds(identifier, nbres, dsspMasking, viterbiparams):
    """ Set the thresholds at which the Viterbi algorithm is stopped """

    if dsspMasking:
        viterbirecalclower = viterbiparams["viterbilowerthreshDSSP"]
        viterbirecalclowertimesaving = viterbiparams["tsviterbilowerthreshDSSP"]
        viterbirecalclowertimesavingplus = viterbiparams["xtsviterbilowerthreshDSSP"]
    else:
        # Save more time here because there are more PSIPRED-based than DSSP-based predictions
        viterbirecalclower = viterbiparams["viterbilowerthreshPSIPRED"]
        viterbirecalclowertimesaving = viterbiparams["tsviterbilowerthreshPSIPRED"]
        viterbirecalclowertimesavingplus = viterbiparams["xtsviterbilowerthreshPSIPRED"]


    if nbres > viterbiparams["tspdbsize"]:
        print("\nWARNING: Lowering the threshold for Viterbi detection because this protein %s is too big (time-saving threshold)"%identifier)
        viterbirecalclower = viterbirecalclowertimesaving
        viterbirecalclowerPSM = viterbirecalclowertimesavingplus
    elif nbres > viterbiparams["xtspdbsize"]:
        print("\nWARNING: Lowering the threshold for Viterbi detection because this protein %s is too big (extra-time-saving threshold)"%identifier)
        viterbirecalclower = viterbirecalclowertimesavingplus
        viterbirecalclowerPSM = viterbirecalclowertimesavingplus
    else:
        viterbirecalclowerPSM = viterbirecalclowertimesaving

    return viterbirecalclower, viterbirecalclowerPSM


def getProbaTransformedGamma(x0, paramshere,paramsbackground, posshift):
    """ Given a coupling x, calculate the log-odds ratio (positive / background) """
    alpha_pos,beta_pos,b_pos,alpha_neg,b_neg,w,shift = [paramshere[p] for p in ["alpha_pos","beta_pos","b_pos","alpha_neg","b_neg","w","shift"]]
    alpha_pos_bg,beta_pos_bg,b_pos_bg,alpha_neg_bg,b_neg_bg,w_bg,shift_bg = [paramsbackground[p] for p in ["alpha_pos","beta_pos","b_pos","alpha_neg","b_neg","w","shift"]]

    x = x0-shift
    xbg = x0-shift_bg
    if x<=0:
        res = np.log((1-w) * b_neg**(-1/alpha_neg)/gamma(1/alpha_neg) * np.fabs(alpha_neg) * np.exp(-np.fabs(x)**alpha_neg/b_neg))
    else:
        res = np.log(w * b_pos**(-1/alpha_pos)/gamma(1/alpha_pos) * np.fabs(alpha_pos/beta_pos) * 1/(1+(x/beta_pos)**alpha_pos)**(1/alpha_pos+1) * np.exp(-(x/beta_pos)**alpha_pos /(b_pos * (1+ (x/beta_pos)**alpha_pos))))
    if xbg <= 0:
        res -= np.log((1-w_bg) * b_neg_bg**(-1/alpha_neg_bg)/gamma(1/alpha_neg_bg) * np.fabs(alpha_neg_bg) * np.exp(-np.fabs(xbg)**alpha_neg_bg/b_neg_bg))
    else:
        res -= np.log(w_bg * b_pos_bg**(-1/alpha_pos_bg)/gamma(1/alpha_pos_bg) * np.fabs(alpha_pos_bg/beta_pos_bg) * 1/(1+(xbg/beta_pos_bg)**alpha_pos_bg)**(1/alpha_pos_bg+1) * np.exp(-(xbg/beta_pos_bg)**alpha_pos_bg /(b_pos_bg * (1+ (xbg/beta_pos_bg)**alpha_pos_bg))))

    # This is to ensure that however the fits for the coupling distributions look like,
    # on the negative side we force the log-odds ratio to be maximum 0 (i.e. we force the background to be above the signal)
    # and in the right-side part of the positive side we force the log-odds ratio to be minimum 0 (i.e. we force the signal to be above the background)
    if xbg < 0 and res > 0:
        res = 0
    if x0 > posshift+0.1 and res < 0:
        res = 0

    return res


def calculateEmission(direction, nbres, predcouplings, secstructseq, secprobdic, params, lowDiversity, i, j, patt):
    """ Calculate the emission probabilities as log-odds ratios for this cell (i,j) """
    dic = {}
    # Because all states except background have the same coupling-based emission probabilities, calculate the product only once
    state="internal"
    posshift = params[direction][state][2]["shift"]

    # calculate the sum of coupling-based log-odds for all cells in the pattern
    probcoupling = 0.0
    if not lowDiversity:
        for m,xy in enumerate(patt):
            if xy=="NA":
                continue
            (x,y) = xy
            coupling = predcouplings[x][y]
            # here we normalized each probability by the background model
            prod = getProbaTransformedGamma(coupling,params[direction][state][m+1],params[direction]["background"][m+1],posshift)
            probcoupling += prod

    if direction=="Antiparallel":
        jprev = j+1
    else:
        jprev = j-1

    liststates = ["bulgei0", "bulgei1", "bulgei2", "bulgej0", "bulgej1", "bulgej2", "first", "internal", "last"]
    for state in liststates:
        # tab is the combination of secondary structure states for this state and this cell (i,j)
        if jprev > 0:
            tab = secstructseq[i-1]+secstructseq[j-1]+secstructseq[i-2]+secstructseq[jprev-1]
        else:
            tab = secstructseq[i-1]+secstructseq[j-1]+secstructseq[i-2]+"C"
        probsec = secprobdic[direction][state][tab]
        # the final emission is the sum of the log-odds parts (secondary-structure-based and coupling-based)
        dic[state] = probsec + (probcoupling)

    return dic


def FindPatternAroundParallel(i, j, size):
    """ Find the cells belonging to the pattern around cell (i,j) in the parallel case """
    tab = []
    for k in range(-1,2):
        # check that we do not run out of the contact map
        if (i+k) > size or (j-k) <= 0:
            tab.append("NA")
            continue
        # check that we do not cross the main diagonal of the contact map
        if (i+k) <= (j-k):
            tab.append("NA")
            continue
        tab.append([i+k,j-k])

    return tab

def FindPatternAroundAntiparallel(i, j, size):
    """ Find the cells belonging to the pattern around cell (i,j) in the antiparallel case """
    tab = []
    for k in range(-1,2):
        # check that we do not run out of the contact map
        if (i+k) > size or (j+k) > size or (i+k) <=0 or (j+k) <=0:
            tab.append("NA")
            continue
        # check that we do not cross the main diagonal of the contact map
        if (i+k) <= (j+k):
            tab.append("NA")
            continue
        tab.append([i+k,j+k])

    return tab


def calculateAllEmissions(nbres, predcouplings, secstructseq, secprobdic, params, lowDiversity):
    """ Calculate all emission probabilities """
    emissions = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(float))))
    patternfunction = {"Parallel":FindPatternAroundParallel, "Antiparallel":FindPatternAroundAntiparallel}
    for direction in ["Parallel","Antiparallel"]:
        for j in range(1,nbres+1):
            for i in range(j+1,nbres+1):
                # for each position, find the cells belonging to the pattern around this cell
                patt = patternfunction[direction](i,j,nbres)
                # calculate the emission probability accordingly
                emissions[direction][i][j] = calculateEmission(direction, nbres, predcouplings, secstructseq, secprobdic, params, lowDiversity, i, j, patt)
    return emissions

