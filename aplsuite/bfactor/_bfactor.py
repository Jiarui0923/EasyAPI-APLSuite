from biopandas.pdb import PandasPdb
import warnings

def get_bfactor(pdb, record='ATOM'):
    warnings.filterwarnings('ignore')
    _pdb_lines = str(pdb).splitlines(True)
    _pdb_df = PandasPdb()._construct_df(_pdb_lines)
    if record == 'ATOM':
        _df = _pdb_df['ATOM'][_pdb_df['ATOM'].atom_name == 'CA']
        _df = _df.groupby(['residue_number']).apply(lambda x: x[x.occupancy==x.occupancy.max()], include_groups=False)
        return _df.b_factor.to_list()
    elif record == 'HETATM':
        _df = _pdb_df['HETATM'][_pdb_df['HETATM'].atom_name == 'CA']
        _df = _df.groupby(['residue_number']).apply(lambda x: x[x.occupancy==x.occupancy.max()], include_groups=False)
        return _df.b_factor.to_list()
    else: raise ValueError(f'{record} Not A Acceptable Record Name.')