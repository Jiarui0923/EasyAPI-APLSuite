from biopandas.pdb import PandasPdb
import warnings

def get_pdb(pdb_id, source='alphafold2-v4'):
    warnings.filterwarnings('ignore')
    if len(pdb_id) == 4: _pdb = PandasPdb().fetch_pdb(pdb_code=pdb_id, source='pdb')
    elif len(pdb_id) == 6: _pdb = PandasPdb().fetch_pdb(uniprot_id=pdb_id, source=source)
    else: raise ValueError('Irregular PDB ID')
    _pdb = _pdb.to_pdb_stream().read()
    return _pdb