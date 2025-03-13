from ..types import Alleles, MHCIIMethods, MHCII, RegularedMers, MHCIIScore, Switcher
from easyapi import register, cache, stat
from ._mhcii import get_mhcii
import json

@register(required_resources={'cpu':1, 'cuda':0})
@stat()
@cache(disable=False)
def mhcii(regular_peptides: RegularedMers['JSON file for peptides and its size and hop.'],
          alleles: Alleles['The alleles for this sequence, seperate by `,`.'] = 'DPA1*01:03',
          mhc_method: MHCIIMethods['The method used to compute MHC-II binding.'] = 'NN_align-NetMHCII-2_3',
          mhc_score: MHCIIScore['The IEDB MHC-II Response Metrics.'] = 'rank',
          fuzzy_alleles: Switcher['Allow alleles to be fuzzy aligned to IEDB format automatically.'] = 'on',
          resources={}) -> dict[
            MHCII['mhc', 'The MHC-II binding outputs from IEDB following JSON formats.'],
            Alleles['alleles', 'The standarized alleles that were sent to IEDB.'],
         ]:
    '''Get MHC-II Binding Prediction
    Use IEDB to predict the MHC-II binding.
    '''
    alleles = [allel for allel in alleles.split(',') if len(allel) > 0]
    fuzzy_alleles = True if fuzzy_alleles == 'on' else False
    _mhcii, _alleles = get_mhcii(mers=list(json.loads(regular_peptides).keys()), method=mhc_method,
                       alleles=alleles, mhc_score=mhc_score, fuzzy_alleles=fuzzy_alleles)
    return dict(mhc=_mhcii, alleles=_alleles)