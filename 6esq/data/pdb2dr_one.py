import sys
import os

import Bio.PDB
import Bio.PDB.StructureBuilder
from Bio.PDB.Residue import Residue
from Bio.PDB import PDBParser

### Distance calculation
parser = PDBParser(PERMISSIVE=1)

# Read structure from files
f = open('pdb2dr_AC.csv','w')

structure_id = "AC"
filename = "crosslinking/6ESQ_AC_A-6ESQ_AC_C.pdb"
structure = parser.get_structure(structure_id, filename)
model = structure[0]

chain_lst = []
for chain in model.get_chains():
    chain_lst.append(chain.get_id())
chainA = model[chain_lst[0]]
chainB = model[chain_lst[1]]
print(chain_lst)

# loop through all residue in chain
# ONLY WORKS FOR DIMER NOW
for residue1 in chainA:
    for residue2 in chainB:
        try:
            # compute distance between CA atoms
            # should use NeighborSearch modeule in biopython
            distance = residue1['CA'] - residue2['CA'] 
        except KeyError:
            ## no CA atom, e.g. for H_NAG
            continue
        if distance<12:
            line = f"6esq_{chainA.get_id()},{residue1.get_id()[1]},6esq_{chainB.get_id()},{residue2.get_id()[1]-592},{distance} \n"
            f.write(line)

f.close()