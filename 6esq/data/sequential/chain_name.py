### only for sequential modeling method

from Bio.PDB import PDBIO, PDBParser

io = PDBIO()
parser = PDBParser()
structure = parser.get_structure('20220721_6esq_aligned', '20220721_6esq_aligned.pdb')

renames = {
    "A": "a",
    "B": "c",
    "C": "h",
    "D": "l",
    "E": "j",
    "F": "e",
    "G": "d",
    "H": "b",
    "I": "g",
    "J": "k",
    "K": "i",
    "L": "f"
}

for model in structure:
    for chain in model:
        old_name = chain.get_id()
        new_name = renames.get(old_name)
        if new_name:
            print(f"renaming chain {old_name} to {new_name}")
            chain.id = new_name
        else:
            print(f"keeping chain name {old_name}")

renames = {
    "a": "A",
    "b": "B",
    "c": "C",
    "d": "D",
    "e": "E",
    "f": "F",
    "g": "G",
    "h": "H",
    "i": "I",
    "j": "J",
    "k": "K",
    "l": "L"
}

for model in structure:
    for chain in model:
        old_name = chain.get_id()
        new_name = renames.get(old_name)
        if new_name:
            print(f"renaming chain {old_name} to {new_name}")
            chain.id = new_name
        else:
            print(f"keeping chain name {old_name}")

io.set_structure(structure)
io.save('20220721_6esq_aligned_renamed.pdb')