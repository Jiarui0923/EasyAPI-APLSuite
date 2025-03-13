from .aplmhc_combine import combine_APL_MHC
import pandas as pd
import io
from pandas.api.types import is_float_dtype
from ..types import PeptideLikelihood, MHCII, PositiveNumber, MHCAPLCombined, Number
from easyapi import register, cache

@register(required_resources={'cpu':1, 'cuda':0})
@cache(disable=True)
def combine_apl_mhc(peptide_likelihood: PeptideLikelihood['Peptide Level Likelihood.'],
                    mhc: MHCII['The MHC-II binding outputs from IEDB following JSON formats.'],
                    w_apl: PositiveNumber['The weight for APL. The weight for MHC will be `1-w_apl`.'] = 0.31,
                    apl_threshold: Number['The threshold for APL. If the value is smaller than this threshold, it will be ignored.'] = 0,
                    mhc_threshold: Number['The threshold for MHC. If the value is smaller than this threshold, it will be ignored.'] = 0,
                    resources = {}) -> dict[
                        MHCAPLCombined['apl_mhc_combined', 'Combined APL-MHC values for each MHC class.']
                    ]:
    '''Weighted Combine APL And MHC
    Combine APL and MHC values use a given weight.
    '''
    _mhc_df = pd.read_json(io.StringIO(mhc))
    _outputs = {}
    for column in _mhc_df.columns:
        if not is_float_dtype(_mhc_df[column]): _outputs[column] = _mhc_df[column].values
        else:
            _combined = combine_APL_MHC(structure=peptide_likelihood, MHC=_mhc_df[column].to_list(), r=w_apl, TProc=apl_threshold, TMHC=mhc_threshold)
            _outputs[column] = _combined
    _outputs = pd.DataFrame(_outputs).to_json(index=False, indent=2)
    return dict(apl_mhc_combined=_outputs)
