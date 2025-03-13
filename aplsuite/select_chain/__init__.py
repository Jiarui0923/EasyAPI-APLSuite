from ..types import PDB, Chain
from easyapi import register, cache, stat
from .select_chain import select_chain as _select_chain

@register(required_resources={'cpu':1, 'cuda':0})
@stat()
@cache(disable=True)
def select_chain(pdb: PDB['The input PDB file.'],
                 chain: Chain['The selected protein chains ID.'] = 'A',
                 resources = {}) -> dict[
                     PDB['pdb', 'The output PDB file that only contains selected chains.']
                 ]:
    """Select Chains from PDB File
    Select destinated chains from the given PDB file.
    """
    _pdb_data = _select_chain(pdb=pdb, chain=chain)
    return dict(pdb=_pdb_data)