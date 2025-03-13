from gpucorex import COREX, Peptide
from gpucorex import sampler as corex_sampler
import torch
import os
import tempfile
from biopandas.pdb import PandasPdb
import numpy as np

def reset_id(input_path='', output_path=None):
    pandas_pdb = PandasPdb().read_pdb(input_path)
    df = pandas_pdb.df['ATOM']
    idx = np.unique(list(zip(df['residue_number'], df['insertion'])),
                    axis=0,
                    return_index=True)[-1]
    idx = sorted(idx)
    _updated_indices = []
    for _i, _o_i in enumerate(idx):
        res_num, res_insert = df.iloc[_o_i].residue_number, df.iloc[_o_i].insertion
        _size = len(df[(df.residue_number==res_num) & (df.insertion==res_insert)])
        _updated_indices +=  ([_i + 1] * _size)
    df.residue_number = np.array(_updated_indices)
    pandas_pdb.df['ATOM'] = df
    
    if output_path is None: output_path = input_path
    pandas_pdb.to_pdb(output_path, records=['ATOM'])
    return pandas_pdb

def corex(pdb:str,
          window_size:int=10,
          min_size:int=4,
          samples:int=10000,
          sampler:str='exhaustive',
          threshold:float=0.75,
          sconf_weight:float=1.0,
          base_fraction:float=1.0,
          probe_radius:float=1.4,
          point_number:float=1000,
          worker_num:int=-1,
          gpu_num:int=-1):
     
     if worker_num <= 0: worker_num=os.cpu_count()-1
     if gpu_num == 0: device='cpu'
     else:
          if gpu_num < 0: device = 'cuda'
          else: device = [f'cuda:{i}' for i in range(0, gpu_num)]
     _pdb_path = ''
     with tempfile.NamedTemporaryFile('w', delete=False, suffix='.pdb') as _pdb_file:
          _pdb_file.write(pdb)
          _pdb_path = _pdb_file.name
     reset_id(_pdb_path)
     _protein = Peptide(path=_pdb_path, window_size=window_size, min_size=min_size)
     try: os.remove(_pdb_path)
     except: pass
     if len(_protein.residues) / window_size < np.log2(samples):
          sampler = 'exhaustive'
     if sampler == 'exhaustive':
          if len(_protein.residues) > 200: raise ValueError('Protein is too large to run exhaustive sampling.')
          sampler = corex_sampler.exhaustive
          sampler_args = {}
     elif sampler == 'montecarlo':
          sampler = corex_sampler.montecarlo
          sampler_args = {'probability': threshold}
     elif sampler == 'adaptive':
          sampler = corex_sampler.adaptive_montecarlo
          sampler_args = {'probability': threshold, 'adaptive_rate':0.05}
     else: raise ModuleNotFoundError(f'Sampler {sampler} Not Found.')
     
     batch_size = 1024
     if len(_protein.residues) < 100: batch_size = 1024
     elif len(_protein.residues) < 200: batch_size = 512
     elif len(_protein.residues) < 400: batch_size = 256
     elif len(_protein.residues) < 600: batch_size = 128
     else: batch_size = 64
     
     while batch_size >= 1:
          try:
               _corex = COREX(workers=worker_num, batch_size=batch_size, samples=samples,
                              device=device, dtype=torch.float64, sampler=sampler, sampler_args=sampler_args,
                              base_fraction=base_fraction, silence=True, probe_radius=probe_radius,
                              point_number=point_number, sconf_weight=sconf_weight)
               return _corex(_protein).tolist()
               break
          except torch.OutOfMemoryError: batch_size = int(batch_size / 2)
          except Exception as e: raise e
     
     