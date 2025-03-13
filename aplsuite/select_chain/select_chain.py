from biopandas.pdb import PandasPdb
import warnings

def _build_chain_condition(df, chains):
    condition = None
    for chain in chains:
        if condition is None: condition = (df.chain_id == chain)
        else: condition = (condition | (df.chain_id == chain))
    return condition

def _build_conditions(df, chains):
    atom_condition = _build_chain_condition(df['ATOM'], chains=chains)
    hetatm_condition = _build_chain_condition(df['HETATM'], chains=chains)
    df['ATOM'] = df['ATOM'][atom_condition]
    df['ATOM'] = df['ATOM'][df['ATOM']['element_symbol'] != 'H']
    df['HETATM'] = df['HETATM'][hetatm_condition]
    df['HETATM'] = df['HETATM'][df['HETATM']['element_symbol'] != 'H']
    df['HETATM'].loc[df['HETATM']['record_name'] == 'HETATM', 'record_name'] = 'ATOM'
    return df

def select_chain(pdb, chain='A'):
    warnings.filterwarnings('ignore')
    _pdb_lines = str(pdb).splitlines(True)
    _pdb_df = PandasPdb()._construct_df(_pdb_lines)
    _chains = set(_pdb_df['ATOM'].chain_id.unique()) | set(_pdb_df['HETATM'].chain_id.unique())
    _target_chains = chain.split(',')
    for _target_chain in _target_chains:
        if _target_chain not in _chains: raise KeyError(f'Chain {_target_chain} Not Found')
    _pdb_df = _build_conditions(_pdb_df, _target_chains)
    _pdb_obj = PandasPdb()
    _pdb_obj._df = _pdb_df
    
    _pdb_io = _pdb_obj.to_pdb_stream(records=['ATOM'])
    return _pdb_io.read()