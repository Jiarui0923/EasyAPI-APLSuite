from ._get_seq import get_seq as get_seq_
from ..types import PDB, Sequence
from easyapi import register, cache

@register(required_resources={'cpu':1, 'cuda':0})
@cache(disable=True)
def get_sequence(pdb: PDB['The input PDB file.'],
                 resources = {}) -> dict[
                     Sequence['sequence', 'The protein amio acid sequence. The order is the same order as the PDB.']
                 ]:
    '''Get Protein Sequence
    Extract the amio acid sequence of the given protein.
    '''
    _seq = get_seq_(pdb=pdb)
    return dict(sequence=_seq)