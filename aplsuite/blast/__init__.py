from ._blast import blast_fasta
import os
from ..types import Sequence, BLASTAlgorithm, BLASTDatabase, BLASTMatrix, FASTA, PositiveNumber, NumberGreaterThan1
from easyapi import register, cache, stat

@register(required_resources={'cpu':-1, 'cuda':0})
@stat()
@cache(disable=False)
def blast(sequence: Sequence['The protein amio acid sequence. The order is the same order as the PDB.'],
          blast_algorithm: BLASTAlgorithm['The BLAST algorithm which could be `blastp` or `blastx`'] = 'blastp',
          blast_database: BLASTDatabase['The BLAST database could be `uniref50`.'] = 'uniref50',
          expect_value: PositiveNumber['The expect threshold sets the maximum e-value threshold for hits to be reported. Lower values make the search more stringent.'] = 10,
          word_size: PositiveNumber['This is the size of initial words or seed matches used in the search. Smaller values increase sensitivity but can slow down the search.'] = 3,
          max_target_seqs: NumberGreaterThan1['Specifies the maximum number of aligned sequences to return. Increasing this will yield more hits.'] = 500,
          matrix: BLASTMatrix['Different matrices can affect the sensitivity for detecting homologous sequences.'] = 'BLOSUM62',
          resources={}) -> dict[
              FASTA['alignment', 'The BLAST outpus in FASTA format.']
          ]:
    '''BLAST Protein Sequence
    BLAST a provided protein sequence and return a FASTA formatted result.
    '''
    _seq = blast_fasta(sequence=sequence, algorithm=blast_algorithm, db=blast_database,
                       num_worker=resources.get('cpu', 1), expect_value=expect_value,
                       word_size=int(word_size), max_target_seqs=int(max_target_seqs),
                       matrix=matrix)
    return dict(alignment=_seq)