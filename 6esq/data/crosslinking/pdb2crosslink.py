import sys

import Bio.PDB
import Bio.PDB.StructureBuilder
from Bio.PDB.Residue import Residue
from Bio.PDB import PDBParser


### Distance calculation
parser = PDBParser(PERMISSIVE=1)

# Read structure from files
structure_id = "6esq_AH"
filename = "6ESQ_AH_A-6ESQ_AH_H.pdb"
structure = parser.get_structure(structure_id, filename)

model = structure[0]
chainA = model['A']
chainB = model['H']

f = open('pdb2cl_AH.csv','w')
pdb2cldb = "Protein 1,Residue 1,Protein 2,Residue 2,UniqueID \n"
f.write(pdb2cldb)

# https://bioinformatics.stackexchange.com/questions/783/how-can-we-find-the-distance-between-all-residues-in-a-pdb-file
# loop through all residue in chain
# ONLY WORKS FOR DIMER NOW
for residue1 in chainA:
    for residue2 in chainB:
        # compute distance between CA atoms
        try:
            distance = residue1['CA'] - residue2['CA']
        except KeyError:
                ## no CA atom, e.g. for H_NAG
                continue
        if distance < 6:
            #print(chainA, residue1, chainC, residue2, distance)
            line = f"6esq_A,{residue1.get_id()[1]},6esq_H,{residue2.get_id()[1]-len(chainA)} \n"
            pdb2cldb += line
            f.write(line)

print(pdb2cldb)
f.close()