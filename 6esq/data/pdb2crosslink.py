from __future__ import print_function

# Biopython
import sys
import Bio.PDB
import Bio.PDB.StructureBuilder
from Bio.PDB.Residue import Residue
from Bio.PDB import PDBParser

# IMP
import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.restraints
import IMP.pmi.restraints.crosslinking
import ihm.cross_linkers

### Distance calculation
parser = PDBParser(PERMISSIVE=1)

# Read structure from files
structure_id = "6esq_AC"
filename = "6ESQ_AC_A-6ESQ_AC_C.pdb"
structure = parser.get_structure(structure_id, filename)

model = structure[0]
chainA = model['A']
chainC = model['C']

cldb = "Protein 1,Protein 2,Residue 1,Residue 2,UniqueID,Score \n"

# loop through all residue in chain
for residue1 in chainA:
    for residue2 in chainC:
        # compute distance between CA atoms
        try:
            distance = residue1['CA'] - residue2['CA']
        except KeyError:
                ## no CA atom, e.g. for H_NAG
                continue
        if distance < 6:
            #print(chainA, residue1, chainC, residue2, distance)
            cldb += f"ProtA,ProtB,{residue1.get_id()[1]},{residue2.get_id()[1]},,1.0 \n"

print(cldb)

### Cross linking data generation
xldb='''Protein 1,Protein 2,Residue 1,Residue 2,UniqueID,Score
ProtA,ProtB,1,10,1,1.0
ProtA,ProtB,1,11,2,2.0
ProtA,ProtB,1,21,3,2.0
'''