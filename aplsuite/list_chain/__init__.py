from ._list_chain import list_chain as list_chain_
from ..types import Chain, PDB
from easyapi import register, cache

@register(required_resources={'cpu':1, 'cuda':0})
@cache(disable=True)
def list_chain(pdb: PDB['The input PDB file.'],
               resources = {}) -> dict[
                   Chain['chain', 'The chains contained in the PDB files.']
                   ]:
    '''List Chains from PDB File
    List all chains from the given PDB file.
    '''
    _chains = list_chain_(pdb=pdb)
    return dict(chain=_chains)