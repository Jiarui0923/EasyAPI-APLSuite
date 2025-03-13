from ._apl import get_apl_peptide
from ..types import Entropy, BFactor, COREX, SASA
from ..types import PeptideLikelihood, ResidueLikelihood, RegularedMers
from ..types import PositiveNumber, NumberGreaterThan1, Number
from easyapi import register, cache, stat
import json

@register(required_resources={'cpu':1, 'cuda':0})
@stat()
@cache(disable=True)
def apl_peptide(residue_likelihood: ResidueLikelihood['Residue Level Likelihood.'],
                regular_peptides: RegularedMers['JSON file for peptides and its size and hop.'],
                resources = {}) -> dict[
                    PeptideLikelihood['peptide_likelihood', 'Peptide Level Likelihood.'],
                ]:
    '''(APL) Peptide-Level Antigen Processing Likelihood
    Calculate peptide level APL based on given residue APL and Peptides.
    '''
    
    peptide_likelihood = get_apl_peptide(residue_likelihood, json.loads(regular_peptides))
    return dict(peptide_likelihood=peptide_likelihood)
