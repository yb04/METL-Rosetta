import requests
import pandas as pd

pdb_id = "1A3A"


def get_AA(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    filePDB = requests.get(url).text.splitlines()
    seq = []
    for line in filePDB:
        if line.startswith("SEQRES"):
            if str(line[11]) == "A":
                seq.append(str(line[19:].strip()))
    result = convert_to_single(seq)
    return result


def convert_to_single(seqeunce):
    bitch = []
    for entry in seqeunce:
        parts = entry.split(" ")
        for part in parts:
            bitch.append(three_to_one_amino_acid_code(part))
    string = "".join(bitch)
    return string


def three_to_one_amino_acid_code(three_letter_code):
    AA_MAPPING = {
        # Standard 20
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLU": "E",
        "GLN": "Q",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",

        # Ambiguity codes (PDB-standard)
        "ASX": "B",  # ASP or ASN
        "GLX": "Z",  # GLU or GLN
        "XLE": "J",  # LEU or ILE
        "XAA": "X",  # unknown generic aa

        # Rare natural amino acids
        "SEC": "U",  # Selenocysteine
        "PYL": "O",  # Pyrrolysine

        # Extremely common PDB modifications
        "MSE": "M",  # Selenomethionine -> treat as Met
        "HYP": "P",  # Hydroxyproline -> Pro
        "PCA": "E",  # Pyroglutamic acid -> Glu
        "CSO": "C",  # S-hydroxycysteine -> Cys
        "MLY": "K",  # N-dimethyl-lysine -> Lys
        "DAL": "A",  # D-Alanine
        "DAR": "R",
        "DSG": "N",  # D-Asn
        "DGN": "Q",  # D-Gln
        "DLY": "K",  # D-Lys
        "DPN": "F",  # D-Phe
        "DPR": "P",  # D-Pro
        "DSN": "S",  # D-Ser
        "DTH": "T",
        "DTR": "W",
        "DTY": "Y",
        "DVA": "V",

        # Less common but still realistic
        "SEP": "S",  # Phosphoserine
        "TPO": "T",  # Phosphothreonine
        "PTR": "Y",  # Phosphotyrosine
        "KCX": "K",  # Carboxylated lysine
        "HIC": "H",  # 4-Methyl-histidine
        "CME": "C",  # Modified cysteine
        "OCS": "C",
        "CSX": "C",

        # Unknown
        "UNK": "X",
    }
    upper_code = three_letter_code.strip().upper()
    one_letter_code = AA_MAPPING.get(upper_code, 'X')  # X für unbekannt
    return one_letter_code


def apply_mutations_ignore(seq, muts):
    if isinstance(muts, str):
        muts = [muts]
    seq_list = list(seq)
    for m in muts:
        pos = int(m[1:-1])
        alt = m[-1]
        seq_list[pos-1] = alt
    return "".join(seq_list)


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

