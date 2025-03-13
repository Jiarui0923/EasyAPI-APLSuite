from biopandas.pdb import PandasPdb
import warnings

def list_chain(pdb):
    warnings.filterwarnings('ignore')
    _pdb_lines = str(pdb).splitlines(True)
    _pdb_df = PandasPdb()._construct_df(_pdb_lines)
    _chains = set(_pdb_df['ATOM'].chain_id.unique()) | set(_pdb_df['HETATM'].chain_id.unique())
    _chains = sorted(_chains)
    _chains = ''.join([f'{_chain},' for _chain in _chains])[:-1]
    return _chains