from scipy import stats
import math
import numpy as np

def smoother(X):
    return (1.0 * sum(X))/len(X)

# reduce window size at ends? currently we do no smoothing in the first and last k-1 residues.
def smooth(L, k):
    if (k % 2 == 0):
        print("Error: window size must be odd!")
        return
    r = int((k-1)/2)
    smoothed_L = L[0:r]
    for i in range(r, len(L) - r):
        window_start = i - r
        window_end = i + r + 1
        smoothed_val = smoother(L[window_start:window_end])
        smoothed_L.append(smoothed_val)
    smoothed_L += L[(len(L)-r):]
    return smoothed_L

#the method combine_zscores is called within the aggregator method
def combine_zscores(COREX, ASA, bfactors, seq_cons, w):
    combined = len(COREX)*[0.0]
    for i in range(len(COREX)):
        combined[i] = (w[0]*seq_cons[i] + w[1]*bfactors[i] + w[2]*COREX[i] + w[3]*ASA[i])/float(math.sqrt(sum(np.square(w))))
    return combined

def filter(x,window_len=10,window='hanning'):

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = x
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')
    return y

def aggregator(seq_cons, bfactors, COREX, ASA, w):
    ASA = [-x for x in ASA]
    bfactors = [-x for x in bfactors]
    seq_cons = [-x for x in seq_cons]

    FILTERING_ON = True
    SMOOTHING_ITERATIONS = 2
    SMOOTHING_WINDOW = 15
    # window size for b-factor?
    for i in range(SMOOTHING_ITERATIONS):
        bfactors = smooth(bfactors, SMOOTHING_WINDOW)
    # window size for sequence entropy? how to get sequence entropy? CLUSTAL-W -> BioEdit
    for i in range(SMOOTHING_ITERATIONS):
        seq_cons = smooth(seq_cons, SMOOTHING_WINDOW)
    # was solvent accessibility smoothed?
    for i in range(SMOOTHING_ITERATIONS):
        ASA = smooth(ASA, SMOOTHING_WINDOW)
    for i in range(SMOOTHING_ITERATIONS):
        COREX = smooth(COREX, SMOOTHING_WINDOW)

    zCOREX = [x for x in stats.zscore(COREX)]
    if (math.isnan(sum(zCOREX))):
        zCOREX = len(zCOREX)*[0.0]
    zASA = [x for x in stats.zscore(ASA)]
    if (math.isnan(sum(zASA))):
        zASA = len(zASA)*[0.0]
    zbfactors = [x for x in stats.zscore(bfactors)]
    if (math.isnan(sum(zbfactors))):
        zbfactors = len(zbfactors)*[0.0]
    zseq_cons = [x for x in stats.zscore(seq_cons)]
    if (math.isnan(sum(zseq_cons))):
        zseq_cons = len(zseq_cons)*[0.0]

    aggregate_zscore = [x for x in combine_zscores(zCOREX, zASA, zbfactors, zseq_cons, w)]

    if(FILTERING_ON == True):
        aggregate_zscore = np.array(aggregate_zscore)
        aggregate_zscore = filter(aggregate_zscore)


    return aggregate_zscore, zCOREX, zASA, zbfactors, zseq_cons
def domain_identifier(raw_scores, LOOP_SIZE = 10, LOWER_THRESHOLD = 0):

    unstable_indices = []
    stable_indices = []
    for i in range(len(raw_scores)):
        if ((raw_scores[i]) >= LOWER_THRESHOLD):
            stable_indices.append(i)
        else:
            unstable_indices.append(i)

    unstable = []
    curr_unstable = [unstable_indices[0], 1, "UNSTABLE"]
    for x in unstable_indices[1:]:
        # the current index is a continuation of the current loop
        if (x == curr_unstable[0]+curr_unstable[1]):
            curr_unstable[1] += 1
        # otherwise we are starting a new loop
        else:
            unstable.append(curr_unstable)
            curr_unstable = [x, 1, "UNSTABLE"]
    unstable.append(curr_unstable)

    stable = []
    curr_stable = [stable_indices[0], 1, "STABLE"]
    for x in stable_indices[1:]:
        # the current index is a continuation of the current domain
        if (x == curr_stable[0]+curr_stable[1]):
            curr_stable[1] += 1
        # otherwise we are starting a new domain
        else:
            stable.append(curr_stable)
            curr_stable = [x, 1, "STABLE"]
    stable.append(curr_stable)

    new_unstable = []
    for l in unstable:
        if (l[1] < LOOP_SIZE):
            stable.append([l[0], l[1], "STABLE"])
        else:
            new_unstable.append(l)

    stable = sorted(stable)
    
    i = 0
    while (i<len(stable)-1):
        if (stable[i][0]+stable[i][1] == stable[i+1][0]):
            stable[i][1] = stable[i][1]+stable[i+1][1]
            stable.pop(i+1)
        else:
            i += 1

    unstable = new_unstable

    X = [max(x, LOWER_THRESHOLD) for x in raw_scores]
    Xnew = X[:]

    return Xnew, stable, unstable


def join_point(segment, i, j):
    i = int(i)
    j = int(j)
    if (i < 0 or j >= len(segment)):
        return
    increment = 1.0*(segment[j]-segment[i])/(j-i)
    k = 0
    while (k+i < j):
        segment[k+i] = segment[i]+k*increment
        k += 1

def weighing_method(Xnew, stable, unstable, MAG , FLANK_SIZE):

    MINIMUM_FLANK_SIZE = 1
    for i in range(len(stable)):
        stable_start = stable[i][0]
        stable_size = stable[i][1]
        stable_end = stable_start+stable_size-1

        if (stable_size < MINIMUM_FLANK_SIZE*2):
            curr_flank_size = min(FLANK_SIZE,int(stable_size*0.5))
            stable[i].append(int(stable_start+(curr_flank_size/2)))
            stable[i].append(int(stable_end-(curr_flank_size/2)))

            for i in range(stable_size):
                Xnew[stable_start+i] = Xnew[stable_start+i]*MAG

        else:
            curr_flank_size = min(FLANK_SIZE,int(stable_size*0.5))

            stable[i].append(int(stable_start+(curr_flank_size/2)))
            stable[i].append(int(stable_end-(curr_flank_size/2)))

            for i in range(curr_flank_size):
                Xnew[stable_start+i] = Xnew[stable_start+i]*MAG
                Xnew[stable_end-i] = Xnew[stable_end-i]*MAG

    segments = sorted(unstable+stable)

    i = 0
    if (len(segments) == 1):
        return Xnew,segments
    while (i < len(segments)):
        if (segments[i][2] == "STABLE"):
            i += 1
        else:
            loop_midpoint = segments[i][0] + int(segments[i][1]/2)

            # interpolate left as long as we have a stable segment to the left
            if ((i > 0) and (segments[i-1][2] == "STABLE")):
                join_point(Xnew, segments[i-1][4], loop_midpoint)
            # interpolate right as long as we have a stable segment to the right
            if ((i < len(segments) - 1) and (segments[i+1][2] == "STABLE")):
                join_point(Xnew, loop_midpoint, segments[i+1][3])
            i += 1

    return Xnew,segments

def residue_APL(seq_cons, bfactors, COREX, ASA,
                w = [1.0,1.0,1.0,1.0], MAG = 2.0, FLANK_SIZE = 24, LOOP_SIZE=10):
    aggregate_zscores,zCOREX, zASA, zbfactors, zseq_cons = aggregator(seq_cons, bfactors, COREX, ASA, w)
    thresholded_zscores, stable, unstable = domain_identifier(aggregate_zscores, LOOP_SIZE)
    residue_scores,segments = weighing_method(thresholded_zscores, stable, unstable, MAG, FLANK_SIZE)
    return [float(i) for i in residue_scores], aggregate_zscores.tolist()

def normalize(list):
    normalized_list = []
    for i in range(len(list)):
        normalized_list.append(float(((list[i]-np.min(list)))/(np.max(list)-np.min(list))))
    return normalized_list

def peptide_APL(residue_scores, mer_size=15, hop=7):

    peptide_likelihood = []
    peptide_pos = []
    for i in range(0, len(residue_scores) - mer_size + 1, hop):
        peptide_pos.append((i, i+mer_size-1))
    for i in range(len(peptide_pos)):
        peptide_likelihood.append(np.mean((residue_scores)[peptide_pos[i][0]:peptide_pos[i][1]+1]))
    peptide_likelihood = normalize(peptide_likelihood)
    return peptide_likelihood

def get_apl_residue(seq_entropy, bfactor, corex, sasa,
            flank_size=20, loop_size=21,
            w_seq_entropy=0.3474973544973545,
            w_bfactor=0.1643121693121693,
            w_corex=0.2651851851851852,
            w_sasa=0.22300529100529098):
    residue_likelihood, aggregate_zscores = residue_APL(seq_entropy, bfactor, corex, sasa,
                                                        FLANK_SIZE=flank_size, LOOP_SIZE=loop_size,
                                                        w=[w_seq_entropy, w_bfactor, w_corex, w_sasa])
    return residue_likelihood, aggregate_zscores
    