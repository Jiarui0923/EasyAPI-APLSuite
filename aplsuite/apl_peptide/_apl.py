import numpy as np

def normalize(list):
    normalized_list = []
    for i in range(len(list)):
        normalized_list.append(((list[i]-min(list)))/(max(list)-min(list)))
    return normalized_list

def peptide_APL(residue_scores, mers):
    peptide_likelihood = []
    pos = 0
    for prop in mers.values():
        pos += (0 if prop.get('step') is None else prop.get('step'))
        _vals = residue_scores[pos:pos+prop['length']]
        peptide_likelihood.append(sum(_vals)/len(_vals))
    peptide_likelihood = normalize(peptide_likelihood)
    return peptide_likelihood

def get_apl_peptide(residue_likelihood, mers):
    peptide_likelihood = peptide_APL(residue_likelihood, mers)
    return peptide_likelihood
    