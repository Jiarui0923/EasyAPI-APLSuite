id='list-chain'
name='List Chains from PDB File'
description='List all chains from the given PDB file.'
version='0.0.1'
references=[]
required_resources={'cpu':1, 'cuda':0}
_types = {
    'pdb': dict(id='pdb',
                meta='string',
                name='PDB File',
                doc='The protein PDB file',
                condition=None,
                version='0.0.1'),
    'chain': dict(id='chain-ids',
                  meta='string',
                  name='PDB Chain IDs',
                  doc='The protein chain ids, seperate with `,`, no blank character.',
                  condition='[A-Za-z0-9]+(,[A-Za-z0-9]+)*',
                  version='0.0.1'),
}
in_params={
    'pdb':   dict(io_type=_types['pdb'], default_value=None, desc='The input PDB file.'),
}
out_params={
    'chain': dict(io_type=_types['chain'], default_value=None, desc='The chains contained in the PDB files.'),
}

from algorithms.proteins.list_chain._list_chain import list_chain
def main(pdb, resources={}):
    _chains = list_chain(pdb=pdb)
    return dict(chain=_chains)
