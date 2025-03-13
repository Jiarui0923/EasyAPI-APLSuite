from biopandas.pdb import PandasPdb
import warnings


def get_sequence(pdb_df):
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    }
    _pdb_res = pdb_df['ATOM'].drop_duplicates(['chain_id','residue_number'])
    _res_seq = _pdb_res.residue_name.apply(lambda x: aa_map.get(x, '?')).values
    _res_seq = ''.join(_res_seq)
    return _res_seq

def get_seq(pdb):
    warnings.filterwarnings('ignore')
    _pdb_lines = str(pdb).splitlines(True)
    _pdb_df = PandasPdb()._construct_df(_pdb_lines)
    return get_sequence(_pdb_df)