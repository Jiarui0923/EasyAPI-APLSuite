import requests
import pandas as pd
import io
import requests
from Levenshtein import distance

species_dict = {
    'Human, HLA-DP': 'http://tools.iedb.org/mhcii/3411/DP/ab/',
    'Human, HLA-DQ': 'http://tools.iedb.org/mhcii/3411/DQ/ab/',
    'Human, HLA-DR': 'http://tools.iedb.org/mhcii/3411/DR/ab/',
    'mouse H-2-I': 'http://tools.iedb.org/mhcii/3411/H2/ab/'
}
def get_alleles(species_dict, species):
    url = species_dict[species]
    content = requests.get(url).content
    decoded_content = content.decode()
    alleles = decoded_content.split("\n")
    return alleles
def make_allele_dict(species_dict):
    alleles_dict = {}
    for species in species_dict:
        alleles_dict[species] = get_alleles(species_dict, species)
    return alleles_dict

def standard_allele():
    alleles_dict = make_allele_dict(species_dict)
    alleles_list = []
    for i in alleles_dict.values(): alleles_list += i
    return alleles_list
    
def allele_translate(alleles):
    regular_alleles = []
    alleles_list = standard_allele()
    for allele in alleles:
        dists = []
        for std_allele in alleles_list:
            dist = distance(allele, std_allele)
            dists.append(dist)
        mapped = alleles_list[dists.index(min(dists))]
        regular_alleles.append(mapped)
    return regular_alleles

def iedb_mhcii(peptides, method='recommended', mhc_score='rank', fuzzy_allele=True,
               alleles=["HLA-DRB1*03:01", "HLA-DRB1*07:01", "HLA-DRB1*15:01", "HLA-DRB3*01:01", "HLA-DRB3*02:02", "HLA-DRB4*01:01", "HLA-DRB5*01:01"]):

    method_mapping = {
        'NN_align-NetMHCII-2_2': 'nn_align-2.2',
        'NN_align-NetMHCII-2_3': 'nn_align-2.3',
        'SMM_align-NetMHCII-1_1': 'smm_align-1.1',
        'NetMHCIIpan-EL-4_1': 'netmhciipan_el-4.1',
        'NetMHCIIpan-EL-4_2': 'netmhciipan_el-4.2',
        'NetMHCIIpan-EL-4_3': 'netmhciipan_el-4.3',
        'NetMHCIIpan-BA-4_1': 'netmhciipan_ba-4.1',
        'NetMHCIIpan-BA-4_2': 'netmhciipan_ba-4.2',
        'NetMHCIIpan-BA-4_3': 'netmhciipan_ba-4.3',
        'NetMHCIIpan-3_2': 'netmhciipan-3.2',
        'NetMHCIIpan-3_1': 'netmhciipan-3.1',
        'Consensus-2_22': 'consensus-2.22',
        'Sturniolo': 'tepitope-1.0',
    }
    score_key = {
        'nn_align-2.2': 'ic50',
        'nn_align-2.3': 'ic50',
        'smm_align-1.1': 'ic50',
        'netmhciipan_el-4.1': 'score',
        'netmhciipan_el-4.2': 'score',
        'netmhciipan_el-4.3': 'score',
        'netmhciipan_ba-4.1': 'score',
        'netmhciipan_ba-4.2': 'score',
        'netmhciipan_ba-4.3': 'score',
        'netmhciipan-3.2': 'ic50',
        'netmhciipan-3.1': 'ic50',
        'consensus-2.22': 'comblib_score',
        'tepitope-1.0': 'score'
    }
    method = method_mapping[method]
    fasta_sequences = "\n".join([f">peptide{i+1}\n{peptides[i]}" for i in range(len(peptides))])
    if fuzzy_allele: alleles = allele_translate(alleles)
    alleles = ','.join(alleles)
    url = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"
    data = {
        "method": method,
        "sequence_text": fasta_sequences,
        "allele": alleles
    }

    response = requests.post(url, data=data)

    if response.status_code != 200: raise ConnectionError(f"{response.text}")

    mhcii_res = pd.read_csv(io.StringIO(response.text), sep='\t')
    if 'peptide' not in mhcii_res.columns:
        raise ValueError("Expected 'peptide' column not found in response.")
    if mhc_score=='score': pivot_df = mhcii_res.pivot(index='peptide', columns='allele', values=score_key[method])
    else: pivot_df = mhcii_res.pivot(index='peptide', columns='allele', values='rank')
    pivot_df = pivot_df.fillna(0)
    pivot_df = pivot_df.reset_index()

    pivot_df.rename(columns={'peptide': 'Peptide'}, inplace=True)

    peptides_df_ = None
    for p in peptides:
        if p in pivot_df['Peptide'].values:
            if peptides_df_ is None: peptides_df_ = pivot_df[pivot_df['Peptide'] == p]
            else: peptides_df_ = pd.concat([peptides_df_, pivot_df[pivot_df['Peptide'] == p]])
    return peptides_df_.reset_index(drop=True), alleles

def get_mhcii(mers=[], method='recommended',
               alleles=["HLA-DRB1*03:01"], mhc_score='rank', fuzzy_alleles=True):
    _mhcii_data, _alleles = iedb_mhcii(method=method, peptides=mers,
                             alleles=alleles, mhc_score=mhc_score, fuzzy_allele=fuzzy_alleles)
    return _mhcii_data.to_json(index=False, indent=2), _alleles