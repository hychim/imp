conda activate ./imp_conda
python script/pdb2dr.py --dimer data/dr_data --dist 6 --outdir data/dr.csv

python script/modeling.py --topology data/topology.txt --fasta data/ --pdbdir data/pdb/ \
--dr data/dr.csv --steps 20000 --outdir result/output/
