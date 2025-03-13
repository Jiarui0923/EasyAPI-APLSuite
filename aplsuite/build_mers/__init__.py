from ._build_mers import build_mers_
from ..types import Sequence, NumberGreaterThan1, Mers
from easyapi import register, cache, stat

@register(required_resources={'cpu':1, 'cuda':0})
@stat()
@cache(disable=True)
def build_mers(sequence: Sequence['The protein amio acid sequence. The order is the same order as the PDB.'],
               mer_size: NumberGreaterThan1['The size of each mer of the sequence.'] = 15,
               hop: NumberGreaterThan1['The size of each hop of mers.'] = 7,
               resources = {}) -> dict[
                    Mers['peptides', 'The protein sequence mers.']
               ]:
    '''Build Peptides from Sequence
    Build peptides with given mer size and hop.
    '''
    _mers = build_mers_(sequence=sequence, mer_size=int(mer_size), hop=int(hop))
    return dict(peptides=_mers)