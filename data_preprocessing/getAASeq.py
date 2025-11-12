import requests
import pandas as pd

pdb_id = "1A3A"

def get_AA(pdbId):
    url = f"https://files.rcsb.org/download/{pdbId}.pdb"
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
    mapping = {
        'ALA': 'A',
        'ARG': 'R',
        'ASN': 'N',
        'ASP': 'D',
        'CYS': 'C',
        'GLU': 'E',
        'GLN': 'Q',
        'GLY': 'G',
        'HIS': 'H',
        'ILE': 'I',
        'LEU': 'L',
        'LYS': 'K',
        'MET': 'M',
        'PHE': 'F',
        'PRO': 'P',
        'SER': 'S',
        'THR': 'T',
        'TRP': 'W',
        'TYR': 'Y',
        'VAL': 'V'
    }
    upper_code = three_letter_code.upper()
    one_letter_code = mapping.get(upper_code, '?')
    return one_letter_code
