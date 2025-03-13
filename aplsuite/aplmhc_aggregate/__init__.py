from ..types import MHCAPLCombined, PeptideLikelihood, MHCII, APLMHCTable, RegularedMers, PositiveNumber, String
from ..types import APLTable
from easyapi import register, cache

import json
import numpy as np

@register(required_resources={'cpu':1, 'cuda':0})
@cache(disable=True)
def apl_mhc_aggregate(peptide_likelihood: PeptideLikelihood['Peptide Level Likelihood.'],
                        apl_mhc_combined: MHCAPLCombined['Combined APL-MHC values for each MHC class.'],
                        mhc: MHCII['The MHC-II binding outputs from IEDB following JSON formats.'],
                        regular_peptides: RegularedMers['JSON file for peptides and its size and hop.'],
                        apl_threshold: PositiveNumber['Threshold for positive.'] = 0.1,
                        antigen: String["Antigen Name Tag"] = '',                
                        resources = {}) -> dict[
                            APLMHCTable['aplmhc_table', 'Combined APL-MHC and its components values.'],
                        ]:
    '''APL-MHC Result Aggregator
    Aggregate APL, MHC, and APL-MHC combined result together to be an intigrated table.
    '''
    data = {
        'MHC': mhc,
        'APL-MHC': apl_mhc_combined,
        'Peptide-Likelihood': peptide_likelihood,
        'Labels': [i['label'] for i in json.loads(regular_peptides).values()],
        'APL-Threshold': apl_threshold,
        'Antigen': antigen
    }
    data = json.dumps(data, indent=2)
    return dict(aplmhc_table=data)
