from python_on_whales import docker
import os
import tempfile
import random
import string
import re

DATABASE_PATH = "/var/services/blast/data/blastdb_custom"

def blast(sequence, algorithm='blastp', db='uniref50',
          num_worker=16, expect_value=10, word_size=3,
          max_target_seqs=500, matrix='BLOSUM62'):
    _path = ''
    with tempfile.NamedTemporaryFile(delete=False) as temp_query:
        temp_id = random.sample(string.ascii_letters, 16)
        temp_query.write(f'>{temp_id}\n{sequence}\n'.encode())
        temp_query.flush()
        
        _path = temp_query.name
    _name = os.path.basename(temp_query.name)
    _query_path_docker = f"/blast/queries/{_name}"
    _data_path_docker = "/blast/blastdb_custom"
    blast_outpus = docker.run(
        "ncbi/blast",
        [algorithm,
            "-query", _query_path_docker,
            "-db", str(db),
            "-num_threads", str(int(num_worker)),
            "-outfmt", "6 sseqid sseq",
            "-evalue", str(expect_value),  # Expect threshold
            "-word_size", str(word_size),  # Word size
            "-max_target_seqs", str(max_target_seqs),  # Max target sequences
            "-matrix", str(matrix),  # Scoring matrix (e.g., BLOSUM62)
            "-max_hsps", "1"  # Max matches in a query range
        ],
        remove=True,
        volumes=[
            (DATABASE_PATH, _data_path_docker, "ro"),
            (_path, _query_path_docker, "ro"),
        ]
    )
    blast_outpus = re.findall(r'(\S+)\s(\S+)', blast_outpus)
    os.remove(_path)
    return blast_outpus
    
def to_fasta(data):
    _fasta = [f'>{_name}\n{_seq}\n' for (_name, _seq) in data]
    return ''.join(_fasta)

def blast_fasta(sequence, algorithm='blastp', db='uniref50',
                num_worker=16, expect_value=10, word_size=3,
                max_target_seqs=500, matrix='BLOSUM62'):
    _blast_out = blast(sequence=sequence, algorithm=algorithm, db=db,
                       num_worker=num_worker, expect_value=expect_value, word_size=word_size,
                       max_target_seqs=max_target_seqs, matrix=matrix)
    return to_fasta(_blast_out)