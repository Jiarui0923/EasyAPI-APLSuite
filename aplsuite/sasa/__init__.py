from ..types import PDB, SASAlgorithm, SASA, NumberGreaterThan1
from easyapi import register, cache, stat
from .sasa_ import sasa as _sasa

@register(required_resources={'cpu':1, 'cuda':0})
@stat()
@cache(disable=False)
def sasa(pdb: PDB['The input PDB file.'],
         algorithm: SASAlgorithm['The SASA algorithms.'] = 'ShrakeRupley',
         probe_radius: NumberGreaterThan1['The probe radius in A.'] = 1.4,
         n_points: NumberGreaterThan1['The number of test points in Shrake & Rupley algorithm.'] = 1000,
         n_slices: NumberGreaterThan1['Get the number of slices per atom in Lee & Richards algorithm.'] = 20,
         resources={}) -> dict[
             SASA['sasa', 'The solvent accessible surface area. The order is the same order as the PDB.']
         ]:
    """Solvent Accessible Surface Area
    Calculate the solvent accessible surface area for the given protein.
    The results will be an array concatenated by the order of sorted(chains).
    """
    sasa_values = _sasa(pdb,
                        algorithm=algorithm,
                        probe_radius=probe_radius,
                        n_points=int(n_points),
                        n_slices=int(n_slices),
                        worker_num=1)
    return dict(sasa=sasa_values)