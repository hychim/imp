import sys
import os
import argparse

import Bio.PDB
import Bio.PDB.StructureBuilder
from Bio.PDB.Residue import Residue
from Bio.PDB import PDBParser

### setting
pdb_parser = PDBParser(PERMISSIVE=1)

parser = argparse.ArgumentParser(description = 'Converting protein interaction pdb file to distant restraint csv.')
parser.add_argument('--dimer', type=str, help = 'Where is the dimer dir')
parser.add_argument('--dist', type=int, help = 'distance for setting distant restraint')
parser.add_argument('--outdir', type=str, help = 'Where to write the result')
args = parser.parse_args()

dimer_dir = args.dimer
dist = args.dist
outdir = args.outdir

# loop through dimer prediction pdb
directory = os.fsencode(dimer_dir)

pdb_lst = []
for file in os.listdir(directory):
     filename = os.fsdecode(file)
     if filename.endswith(".pdb"):
         pdb_lst.append(filename)

# Read structure from files
f = open(outdir,'w')

for pdb in pdb_lst:
    structure_id = pdb
    filename = f"{dimer_dir}/{pdb}"
    structure = pdb_parser.get_structure(structure_id, filename)
    model = structure[0]

    chain_lst = []
    for chain in model.get_chains():
        chain_lst.append(chain.get_id())
    chainA = model[chain_lst[0]]
    chainB = model[chain_lst[1]]
    print(chain_lst)

    for residue in chainB:
        first = residue.get_id()[1]
        break

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
            if distance < dist:
                line = f"chain_{chainA.get_id()},{residue1.get_id()[1]},chain_{chainB.get_id()},{residue2.get_id()[1]-first+1},{distance} \n"
                f.write(line)

f.close()