from ._corex import corex as get_corex
from ..types import PDB, NumberGreaterThan1, COREX, COREXSampler, PositiveNumber
from easyapi import register, cache, stat

@register(required_resources={'cpu':-1, 'cuda':-1})
@stat()
@cache(disable=False)
def corex(pdb: PDB['The input PDB file.'],
          window_size: NumberGreaterThan1['The protein folding unit size. Also, the number of partition schemes.'] = 10,
          min_size: NumberGreaterThan1['The minumum protein folding unit size.'] = 4,
          samples: NumberGreaterThan1['(Ignore for exhaustive sampling) The sample number for each partition scheme. Total sample number=samples*window_size.'] = 10000,
          sampler: COREXSampler['The COREX states sampler.'] = 'exhaustive',
          threshold: PositiveNumber['(Ignore for exhaustive sampling) The threshold for the sampler.'] = 0.75,
          sconf_weight: PositiveNumber['Entropy factor.'] = 1.0,
          base_fraction: PositiveNumber['The base fraction used to sum all COREX (ln_kf) values.'] = 1.0,
          probe_radius: NumberGreaterThan1['The probe radius for SASA in A.'] = 1.4,
          n_points: NumberGreaterThan1['The number of test points in Shrake & Rupley algorithm for SASA.'] = 1000,
          resources = {}) -> dict[
              COREX['corex', 'The COREX values. The order is the same order as the PDB.']
          ]:
    '''(COREX) CORrelation with hydrogen EXchange protection factors
    An algorithm designed to compute comformational stability of a protein. The results will be an array concatenated by the order of sorted(chains)
    '''
    corex_values = get_corex(pdb, window_size=int(window_size), min_size=int(min_size),
                         samples=int(samples), sampler=sampler, threshold=threshold,
                         sconf_weight=sconf_weight, base_fraction=base_fraction,
                         probe_radius=probe_radius, point_number=int(n_points),
                         worker_num=resources.get('cpu', -1), gpu_num=resources.get('cuda', -1))
    return dict(corex=corex_values)