import freesasa
from io import StringIO
from Bio.PDB import PDBParser
from uuid import uuid4

# ===========================================================
# Author: Jai Bansal
# Created Date: 2024-10-20
# Update: Jai Bansal
# Updated Date: 2024-10-23
# Updated Comment: Add customizable paramters
# Description: User can now change algorithm type, probe radius, n_points, n_slices, and n_threads.
# ===========================================================
def sasa(pdb:str,
         algorithm:str='ShrakeRupley',
         probe_radius:float=1.4,
         n_points:int=1000,
         n_slices:int=20,
         worker_num:int=1):
    """
    Run a computation to calculate the solvent-accessible surface area (SASA) for a molecular structure from a PDB file

    Parameters
    ----------
        pdb_path: str
            Path to the PDB file containing the molecular structure to be analyzed.

        algorithm: str, [optional, default: 'Shrake-Rupley']
            The algorithm to use for calculating the solvent-accessible surface area. Common choices are 'Shrake-Rupley' and 'Lee-Richards'.

        probe_radius: float, [optional, default: 1.4]
            The radius of the solvent probe, usually representing the size of a water molecule. Default is 1.4 Å.

        n_points: int, [optional, default: 100]
            Number of points used for sampling the surface. Higher values improve accuracy but increase computational cost.

        n_slices: int, [optional, default: 20]
            Number of slices used for sampling in cylindrical coordinates. This parameter affects the resolution of the SASA calculation.

        n_threads: int, [optional, default: 1]
            Number of threads to use for parallel computation. Increasing this value can speed up calculations by using multiple cores.
        
    Returns
    -------
    dataframe df
        The function returns a pandas DataFrame containing the unique residue IDs and their relative solvent-accessible surface areas (SASA) after processing a PDB structure.
        
    """
    pdb_io = StringIO(pdb)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(str(uuid4()), pdb_io)
    structure = freesasa.structureFromBioPDB(structure)
    sasa_params = freesasa.Parameters({
        'algorithm'    : algorithm, 
        'probe-radius' : probe_radius, 
        'n-points'     : n_points,
        'n-slices'     : n_slices, 
        'n-threads'    : worker_num,
    })
    sasa_values = freesasa.calc(structure, sasa_params).residueAreas()
    sasa_data = {}
    for chain_id, chain_data in sasa_values.items():
        res_data = [item.relativeTotal for item in chain_data.values()]
        sasa_data[chain_id] = res_data
    sasa_values_list = []
    for key in sorted(sasa_data): sasa_values_list += sasa_data[key]
    return sasa_values_list