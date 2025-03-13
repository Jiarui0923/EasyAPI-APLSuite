from ._plddt import get_plddt
from ..types import PDB, PDBRecord, BFactor
from easyapi import register, cache, stat

@register(required_resources={'cpu':1, 'cuda':0})
@stat()
@cache(disable=True)
def plddt(pdb: PDB['The input AlphaFold predicted PDB file.'],
            record: PDBRecord['The AlphaFold PDB record for pLDDT extraction.'] = 'ATOM',
            resources = {}) -> dict[
                BFactor['bfactor', 'The 100-pLDDT. The order is the same order as the PDB.']
            ]:
    '''Extract 100-pLDDT
    Extract residue level 100-pLDDT from the given AlphaFold predicted PDB file (The pLDDT of CA atom).
    '''
    _bfactor = get_plddt(pdb=pdb, record=record)
    return dict(bfactor=_bfactor)
