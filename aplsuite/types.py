from easyapi.annotations import String
from easyapi.annotations import NumArray
from easyapi.annotations import Number
from easyapi.annotations import NumberGreaterThan1
from easyapi.annotations import PositiveNumber
from easyapi.annotations import JSONString

class Sequence(String):
    id        = 'protein-seq'
    name      = 'Protein Amio Acid Sequence'
    doc       = 'The protein amio acid sequence'
    condition = '[ACDEFGHIKLMNPQRSTVWY?]+'

class Mers(String):
    id        = 'protein-peptides'
    name      = 'Protein Sequence Sliced Peptides'
    doc       = 'The protein amio acid sequence sliced as peptides. Each line would be a peptide.'
    condition = r'([ACDEFGHIKLMNPQRSTVWY?]*,?(0|1|true|false|True|False)?\s*)*'
    
class RegularedMers(JSONString):
    id        = 'protein-regulared-peptide'
    name      = 'Protein Sequence Regulared Sliced Peptides'
    doc       = 'The JSON with peptides as keys and hop, mersize, overlap, and label as values.'
class FASTA(String):
    id        = 'fasta'
    name      = 'FASTA Sequence'
    doc       = 'The FASTA sequence file'
    
class Entropy(NumArray):
    id        = 'sequence-entropy'
    name      = 'Sequence Entropy Values'
    doc       = 'Sequence Entropy Values in Sorted Chain ID Order'

class PDB(String):
    id        = 'pdb'
    name      = 'PDB File'
    doc       = 'The protein PDB file.'

class MMCIF(String):
    id        = 'mmcif'
    name      = 'mmCIF File'
    doc       = 'The Macromolecular Crystallographic Information File.'

class PDBID(String):
    id        = 'pdbid'
    name      = 'PDB ID'
    doc       = 'The PDB ID, which can be 4 chars for RCSB or 6 chars for UniProt.'
    condition = '[A-Za-z0-9]{4}|[A-Za-z0-9]{6}'
    
class PDBSource(String):
    id        = 'pdb_source'
    name      = 'PDB Source'
    doc       = '(alphafold2-v3|alphafold2-v4) The PDB fetch source which could be alphafold2-v3 or alphafold2-v4'
    condition = 'alphafold2-v3|alphafold2-v4'

class PDBRecord(String):
    id        = 'pdb-record'
    name      = 'PDB record'
    doc       = '(ATOM|HETATM) PDB record names which could be ATOM or HETATM.'
    condition = '(ATOM|HETATM)'
class Chain(String):
    id        = 'chain-ids'
    name      = 'PDB Chain IDs'
    doc       = 'The protein chain ids, seperate with `,`, no blank character.'
    condition = '[A-Za-z0-9]+(,[A-Za-z0-9]+)*'
    
class SASA(NumArray):
    id        = 'sasa'
    name      = 'SASA Values'
    doc       = 'SASA Values in Sorted Chain ID Order.'
    
class SASAlgorithm(String):
    id        = 'sasa-algorithm'
    name      = 'SASA Algorithm'
    doc       = '(ShrakeRupley|LeeRichards) The SASA Algorithm that could be ShrakeRupley or LeeRichards.'
    condition = '(ShrakeRupley|LeeRichards)'
    
class Alleles(String):
    id        = 'alleles'
    name      = 'The Alleles Marks'
    doc       = 'The alleles marks seperate by `,`. The avaliable options are from iedb.org.'
    version   = '0.0.1'
    
class MHCIIMethods(String):
    id        = 'mhcii-methods'
    name      = 'The IEDB MHC-II Methods'
    doc       = 'The IEDB MHC-II prediction methods.'
    condition = 'NN_align-NetMHCII-2_2|NN_align-NetMHCII-2_3|SMM_align-NetMHCII-1_1|NetMHCIIpan-EL-4_1|NetMHCIIpan-EL-4_2|NetMHCIIpan-EL-4_3|NetMHCIIpan-BA-4_1|NetMHCIIpan-BA-4_2|NetMHCIIpan-BA-4_3|NetMHCIIpan-3_2|NetMHCIIpan-3_1|Consensus-2_22|Sturniolo'
    version   = '0.0.1'
    
class MHCIIScore(String):
    id        = 'mhcii-score'
    name      = 'The IEDB MHC-II Response Metrics'
    doc       = 'The IEDB MHC-II response metrics could be rank or score.'
    condition = '(rank|score)'
    version   = '0.0.1'
    
class MHCII(JSONString):
    id        = 'mhcii'
    name      = 'IEDB MHCII Predictions'
    doc       = 'A JSON file for MHCII prediction results.'
    version   = '0.0.1'
    
class COREXSampler(String):
    id        = 'corex-sampler'
    name      = 'COREX Sampler'
    doc       = '(exhaustive|montecarlo|adaptive) The COREX micro-states sampler, which could be exhaustive enumerate, Monte Carlo, or Adaptibe Monte Carlo sampler.'
    condition = '(exhaustive|montecarlo|adaptive)'
    version   = '0.0.1'
    
class COREX(NumArray):
    id        = 'corex'
    name      = 'COREX (ln(kf)) Values'
    doc       = 'COREX Values in Sorted Chain ID Order'
    version   = '0.0.1'
    
class BLASTAlgorithm(String):
    id        = 'blast-algorithm'
    name      = 'BLAST Algorithm'
    doc       = 'The BLAST algorithm which could be `blastp` or `blastx`'
    condition = '(blastp|blastx)'
    version   = '0.0.1'
    
class BLASTDatabase(String):
    id        = 'blast-database'
    name      = 'BLAST Databases'
    doc       = 'The BLAST database could be `uniref50`.'
    condition = '(uniref50)'
    version   = '0.0.1'
    
class BLASTMatrix(String):
    id        = 'blast-matrix'
    name      = 'BLAST Matrix'
    doc       = 'Different matrices can affect the sensitivity for detecting homologous sequences.'
    condition = '(BLOSUM45|BLOSUM50|BLOSUM62|BLOSUM80|BLOSUM90|PAM30|PAM70|PAM250)'
    version   = '0.0.1'
    
class BFactor(NumArray):
    id        = 'bfactor'
    name      = 'B-Factor Values'
    doc       = 'B-Factor values in given PDB file atom orders.'
    version   = '0.0.1'
    
class PeptideLikelihood(NumArray):
    id        = 'apl-peptide-likelihood'
    name      = 'Peptide Antigen Processing Likelihood'
    doc       = 'Peptide Level Antigen Processing Likelihood.'
    version   = '0.0.1'

class ResidueLikelihood(NumArray):
    id        = 'apl-residue-likelihood'
    name      = 'Residue Antigen Processing Likelihood'
    doc       = 'Residue Level Antigen Processing Likelihood.'
    version   = '0.0.1'

class APLAggregate(NumArray):
    id        = 'apl-aggregate'
    name      = 'Residue Level Aggregated Score'
    doc       = 'Residue level weight-aggregated COREX, SASA, B-Factor, and Sequence Entropy score.'
    version   = '0.0.1'

class MHCAPLCombined(JSONString):
    id        = 'apl-mhc-combined'
    name      = 'Combined APL-MHC Values'
    doc       = 'Combined APL-MHC values for each MHC class.'
    version   = '0.0.1'
    
class APLTable(JSONString):
    id        = 'apl-table'
    name      = 'Aggregated APL Table'
    doc       = 'Combined APL and its components values.'
    version   = '0.0.1'
    
class APLMHCTable(JSONString):
    id        = 'aplmhc-table'
    name      = 'Aggregated APL-MHC Table'
    doc       = 'Combined APL-MHC and its components values.'
    version   = '0.0.1'
    
class Switcher(String):
    id        = 'switcher'
    name      = 'Switcher Variable'
    doc       = 'Turn on or turn off a variable.'
    condition = '(on|off)'
    version   = '0.0.1'