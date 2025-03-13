from ..types import PDB, MMCIF
from easyapi import register, cache, stat

from Bio.PDB import MMCIFParser, PDBIO
from uuid import uuid4
import io
def mmcif_to_pdb_(mmcif):
    parser = MMCIFParser(QUIET=True)
    with io.StringIO(mmcif) as mmcif_io:
        structure = parser.get_structure(str(uuid4()), mmcif_io)
        pdb = PDBIO()
        pdb.set_structure(structure)
        with io.StringIO() as pdb_io:
            pdb.save(pdb_io)
            pdb_io.seek(0)
            return pdb_io.read()

@register(required_resources={'cpu':1, 'cuda':0})
@stat()
@cache(disable=False)
def mmcif_to_pdb(mmcif: MMCIF['The mmCIF format protein file.'],
            resources = {}) -> dict[
                PDB['pdb', 'The fetched PDB file.']
            ]:
    '''mmCIF to PDB
    Convert mmCIF file to PDB file.
    '''
    _pdb = mmcif_to_pdb_(mmcif=mmcif)
    return dict(pdb=_pdb)