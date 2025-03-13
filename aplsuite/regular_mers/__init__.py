from ._regular_mers import measure_frags
from ..types import Mers, RegularedMers
from easyapi import register, cache, stat

import json

@register(required_resources={'cpu':1, 'cuda':0})
@stat()
@cache(disable=True)
def regular_mers(peptides: Mers['The protein sequence peptides.'],
                 resources = {}) -> dict[
                    RegularedMers['regular_peptides', 'JSON file for peptides and its size, hop, and label.']
                 ]:
    '''Regular Given Peptides
    Build peptides with given mer size and hop.
    '''
    frags = [frag.split(',') for frag in peptides.split() if len(frag) > 0]
    labels = [0 if len(frag) <= 1 else frag[1] for frag in frags]
    frags = [frag[0] for frag in frags]
    _mers = json.dumps(measure_frags(frags=frags, labels=labels), indent=2)
    return dict(regular_peptides=_mers)