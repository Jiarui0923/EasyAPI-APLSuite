from ._apl import get_apl_residue as get_apl_residue_
from ..types import Entropy, BFactor, COREX, SASA
from ..types import ResidueLikelihood, APLAggregate
from ..types import NumberGreaterThan1, Number
from easyapi import register, cache, stat

@register(required_resources={'cpu':1, 'cuda':0})
@stat()
@cache(disable=True)
def apl_residue(sequence_entropy: Entropy['The sequence entropy of the given sequence based on the alignments.'],
        bfactor: BFactor['The B-Factor. The order is the same order as the PDB.'],
        corex: COREX['The COREX values. The order is the same order as the PDB.'],
        sasa: SASA['The solvent accessible surface area.'],
        flank_size: NumberGreaterThan1['The flank size of APL.'] = 11,
        loop_size: NumberGreaterThan1['The loop size of APL.'] = 16,
        w_entropy: Number['The weight for entropy.'] = 0.17,
        w_bfactor: Number['The weight for B-factor.'] = 0.18,
        w_corex: Number['The weight for COREX.'] = 0.51,
        w_sasa: Number['The weight for SASA.'] =0.15,
        resources = {}) -> dict[
            ResidueLikelihood['residue_likelihood', 'Residue Level Likelihood.'],
            APLAggregate['apl_aggregate', 'Residue Level Aggregated Score.']
        ]:
    '''(APL) Residue Level Antigen Processing Likelihood
    Calculate residue level APL by Sequence Entropy, B-Factor, SASA, and COREX.
    '''
    residue_likelihood, aggregate_zscores = get_apl_residue_(seq_entropy=sequence_entropy, bfactor=bfactor, corex=corex, sasa=sasa,
                                                                        flank_size=int(flank_size), loop_size=int(loop_size),
                                                                        w_seq_entropy=w_entropy, w_bfactor=w_bfactor, w_sasa=w_sasa, w_corex=w_corex)
    return dict(residue_likelihood=residue_likelihood, apl_aggregate=aggregate_zscores)
