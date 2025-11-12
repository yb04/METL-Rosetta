import requests

pdb_id = "1h2e"

url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
filePDB = requests.get(url).text.splitlines()
seq = []
for line in filePDB:
    if line.startswith("SEQRES"):
        if str(line[11])== "A":
            seq.append(str(line[19:].strip()))

def convert_to_single(seqeunce):
    result = []
    for entry in seqeunce:
        parts = entry.split(" ")
        for part in parts:
            result.append(three_to_one_amino_acid_code(part))
    string = "".join(result)
    return string

def three_to_one_amino_acid_code(three_letter_code):
    mapping = {
        'ALA': 'A',  # Alanin
        'ARG': 'R',  # Arginin
        'ASN': 'N',  # Asparagin
        'ASP': 'D',  # Asparaginsäure
        'CYS': 'C',  # Cystein
        'GLU': 'E',  # Glutaminsäure
        'GLN': 'Q',  # Glutamin
        'GLY': 'G',  # Glycin
        'HIS': 'H',  # Histidin
        'ILE': 'I',  # Isoleucin
        'LEU': 'L',  # Leucin
        'LYS': 'K',  # Lysin
        'MET': 'M',  # Methionin
        'PHE': 'F',  # Phenylalanin
        'PRO': 'P',  # Prolin
        'SER': 'S',  # Serin
        'THR': 'T',  # Threonin
        'TRP': 'W',  # Tryptophan
        'TYR': 'Y',  # Tyrosin
        'VAL': 'V'  # Valin
    }
    upper_code = three_letter_code.upper()
    one_letter_code = mapping.get(upper_code, '?')
    return one_letter_code

result = convert_to_single(seq)
print(result)

def apply_mutations(seq, muts):
    # seq: str mit 1-basiertem Index (biologische Konvention)
    # muts: z.B. "A123V" oder Liste ["K31G","I102W"]
    if isinstance(muts, str):
        muts = [muts]
    seq_list = list(seq)
    for m in muts:
        ref = m[0]
        pos = int(m[1:-1])
        alt = m[-1]
        if seq_list[pos-1] != ref:
            raise ValueError(f"Ref-Mismatch an {pos}: erwartet {ref}, gefunden {seq_list[pos-1]}")
        seq_list[pos-1] = alt
    return "".join(seq_list)


# Beispiel
mut = apply_mutations(result, ["I147H","E184Y","D197M"])
print(mut)
