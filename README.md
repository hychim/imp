# Assembling protein complex with IMP
This project aims to model large protein complex from combination of protein subunits using IMP.

# Computational Requirement
The only requirement is to clone this github repository and set up a conda environment before running the script.
```bash
clone https://github.com/hychim/imp.git
conda env create --prefix ./imp_conda -f environment.yml
```

## Procedure
You can directly run the imp.sh using the 6esq example. If you want to model your own protein complex, please follow below procedures and prepare all the data needed.

1. Put all your dimer pdb files into the data/dr_data directory
2. Put all your fasta file(s) into the data directory
3. Put all your subunit pdb files into the data/pdb directory
4. Create the topology.txt file for the IMP and put it into the data directory.
5. Run imp.sh
