from ..types import Entropy, BFactor, COREX, SASA
from ..types import ResidueLikelihood, APLAggregate, PeptideLikelihood, PositiveNumber, RegularedMers, String
from ..types import APLTable
from easyapi import register, cache

import json
import numpy as np

def zscore(x):
    x = np.array(x)
    z = (x - np.mean(x)) / np.std(x)
    return z.tolist()

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

def regular(x):
    SMOOTHING_WINDOW = 15
    x = smooth(x, SMOOTHING_WINDOW)
    x = zscore(x)
    return x

@register(required_resources={'cpu':1, 'cuda':0})
@cache(disable=True)
def apl_aggregate(sequence_entropy: Entropy['The sequence entropy of the given sequence based on the alignments.'],
        bfactor: BFactor['The B-Factor. The order is the same order as the PDB.'],
        corex: COREX['The COREX values. The order is the same order as the PDB.'],
        sasa: SASA['The solvent accessible surface area.'],
        apl_aggregate: APLAggregate['The size of each mer of the sequence.'],
        residue_likelihood: ResidueLikelihood['The size of each hop of mers.'],
        peptide_likelihood: PeptideLikelihood['Peptide Level Likelihood.'],
        regular_peptides: RegularedMers['JSON file for peptides and its size and hop.'],
        apl_threshold: PositiveNumber['Threshold for positive.'] = 0.1,
        antigen: String["Antigen Name Tag"] = '',
        resources = {}) -> dict[
            APLTable['apl_table', 'Combined APL and its components values.'],
        ]:
    '''APL Result Aggregator
    Aggregate BFactor, SASA, COREX, Sequenc Entropy, and APL together to be an intigrated table.
    '''
    data = {
        'B-Factor': regular(bfactor),
        'SASA': regular(sasa),
        'COREX': regular(corex),
        'Sequence Entropy': regular(sequence_entropy),
        'Aggregate': regular(apl_aggregate),
        'APL': residue_likelihood,
        'Peptide-Likelihood': peptide_likelihood,
        'APL-Threshold': apl_threshold,
        'Regular-Mers': regular_peptides,
        'Antigen': antigen
    }
    data = json.dumps(data, indent=2)
    return dict(apl_table=data)
