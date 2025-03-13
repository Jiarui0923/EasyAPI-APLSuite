import math
import numpy as np
import scipy.stats

def normalize(list):
    normalized_list = []
    for i in range(len(list)):
        normalized_list.append(float(((list[i]-np.min(list)))/(np.max(list)-np.min(list))))
    return normalized_list

def combine_APL_MHC(structure, MHC ,r=0.5, TProc = 0, TMHC = 0):
    combined = []
    if(MHC == None) or (structure == None):
        return combined
    else:
        w = [r, 1-r]
        z_structure = scipy.stats.zscore(structure)
        z_MHC = scipy.stats.zscore(MHC)
        for i in range(len(structure)):
            combined.append((z_structure[i]*w[0] + z_MHC[i]*w[1])/math.sqrt(sum(np.square(w))))

        combined = normalize(combined)
        MHC = normalize(MHC)
        structure = normalize(structure)

        for i in range(len(combined)):
            if structure[i]<TProc or MHC[i]<TMHC:
                combined[i] = 0
        
        return normalize(combined)
