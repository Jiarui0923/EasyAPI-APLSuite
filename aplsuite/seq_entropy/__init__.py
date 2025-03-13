from ..types import Sequence, FASTA, Entropy
from easyapi import register, cache, stat
from .seq_entropy import get_entropy

@register(required_resources={'cpu':1, 'cuda':0})
@stat()
@cache(disable=False)
def sequence_entropy(sequence: Sequence['The protein amio acid sequence. The order is the same order as the PDB.'],
                alignment: FASTA['The BLAST outpus in FASTA format for the given sequence.'],
                resources = {}) -> dict[
                    Entropy['sequence_entropy', 'The sequence entropy of the given sequence based on the alignments.']
                ]:
    """Sequence Entropy
    Compute sequence entropy for the given sequence and BLAST outputs.
    """
    _entropy = get_entropy(sequence=sequence, alignment=alignment)
    return dict(sequence_entropy=_entropy)