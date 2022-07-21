import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.topology
import ihm.cross_linkers

import sys
import csv

# ---------------------------
# Define Input Files
# ---------------------------
datadirectory = "../../data/seq/"
topology_file = datadirectory+"topology_ACHLJEDBGKI.txt"
#topology_file = datadirectory+"topology_ACH.txt"

dr_csv = datadirectory+'pdb2dr_ACHLJEDBGKI.csv'
#dr_csv = datadirectory+'pdb2dr_AC.csv'

# Initialize model
m = IMP.Model()

# Read in the topology file.
# Specify the directory wheere the PDB files, fasta files and GMM files are
topology = IMP.pmi.topology.TopologyReader(topology_file,
                                           pdb_dir=datadirectory,
                                           fasta_dir=datadirectory,
                                           gmm_dir=datadirectory)

# Use the BuildSystem macro to build states from the topology file
bs = IMP.pmi.macros.BuildSystem(m)

# Each state can be specified by a topology file.
bs.add_state(topology)

# Build the system representation and degrees of freedom
root_hier, dof = bs.execute_macro(max_rb_trans=4.0,
                                  max_rb_rot=0.3,
                                  max_bead_trans=4.0,
                                  max_srb_trans=4.0,
                                  max_srb_rot=0.3)

# Randomize the initial configuration before sampling, of only the molecules
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation=50,
                                    verbose=False,
                                    cutoff=5.0,
                                    niterations=100)

outputobjects = []  # reporter objects (for stat files)

# -----------------------------------
# Define Scoring Function Components
# -----------------------------------

# Connectivity keeps things connected along the backbone
mols = IMP.pmi.tools.get_molecules(root_hier)
for mol in mols:
    molname = mol.get_name()
    IMP.pmi.tools.display_bonds(mol)
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
        mol, scale=2.0)
    cr.add_to_model()
    cr.set_label(molname)
    outputobjects.append(cr)

# Excluded Volume Restraint
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         included_objects=root_hier,
                                         resolution=1000)
ev.add_to_model()
outputobjects.append(ev)

# Distant Constraint
with open(dr_csv, newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

for l in data:
    try:
        dr = IMP.pmi.restraints.basic.DistanceRestraint(root_hier,
                                                        (int(l[1]),int(l[1]),l[0],),
                                                        (int(l[3]),int(l[3]),l[2],),
                                                        distancemin = float(l[4])-1,
                                                        distancemax = float(l[4])+1)		
        dr.add_to_model()
    except:
        print(l)
        pass

# --------------------------
# Monte-Carlo Sampling
# --------------------------

# --------------------------
# Set MC Sampling Parameters
# --------------------------
num_frames = 40000
if '--test' in sys.argv:
    num_frames = 100
num_mc_steps = 10

# This object defines all components to be sampled as well as the
# sampling protocol
mc1 = IMP.pmi.macros.ReplicaExchange0(
    m, root_hier=root_hier, monte_carlo_sample_objects=dof.get_movers(),
    output_objects=outputobjects, monte_carlo_temperature=1.0,
    simulated_annealing=True, simulated_annealing_minimum_temperature=1.0,
    simulated_annealing_maximum_temperature=2.5,
    simulated_annealing_minimum_temperature_nframes=200,
    simulated_annealing_maximum_temperature_nframes=20,
    replica_exchange_minimum_temperature=1.0,
    replica_exchange_maximum_temperature=2.5,
    number_of_best_scoring_models=10, monte_carlo_steps=num_mc_steps,
    number_of_frames=num_frames, global_output_directory="output")

# Start Sampling
mc1.execute_macro()

