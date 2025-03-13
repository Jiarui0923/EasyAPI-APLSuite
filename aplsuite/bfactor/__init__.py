from ._bfactor import get_bfactor
from ..types import PDB, PDBRecord, BFactor
from easyapi import register, cache, stat

@register(required_resources={'cpu':1, 'cuda':0})
@stat()
@cache(disable=True)
def bfactor(pdb: PDB['The input PDB file.'],
            record: PDBRecord['The PDB record for B-Factor extraction.'] = 'ATOM',
            resources = {}) -> dict[
                BFactor['bfactor', 'The B-Factor. The order is the same order as the PDB.']
            ]:
    '''Extract B-Factor
    Extract residue level B-Factor from the given PDB file (The B-Factor of CA atom).
    '''
    _bfactor = get_bfactor(pdb=pdb, record=record)
    return dict(bfactor=_bfactor)
