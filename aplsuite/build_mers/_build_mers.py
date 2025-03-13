def build_mers_(sequence, mer_size=15, hop=7):
    peptides = []
    for i in range(0, len(sequence) - mer_size + 1, hop):
        peptide = sequence[i:i + mer_size]
        peptides.append(f'{peptide},0')
    return '\n'.join(peptides)