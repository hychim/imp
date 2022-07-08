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

# ---------------------------
# Define Input Files
# ---------------------------
datadirectory = "../../data/"
topology_file = datadirectory+"topology_all_cl.txt"
#topology_file = datadirectory+"topology_AC_cl.txt"

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

#["A_1","B_1","C_1","D_1","E_1","F_1","G_1","H_1","I_1","J_1","K_1","L_1"]

# Fix all rigid bodies
# First select and gather all particles to fix.
fixed_particles = []
for prot in []:                            # ADD FIXED PROTEIN HERE
    fixed_particles += IMP.atom.Selection(
        root_hier, molecule=prot).get_selected_particles()

# Fix the Corresponding Rigid movers and Super Rigid Body movers using dof
# The flexible beads will still be flexible (fixed_beads is an empty list)!
fixed_beads, fixed_rbs = dof.disable_movers(fixed_particles,
                                            [IMP.core.RigidBodyMover,
                                             IMP.pmi.TransformMover])

# Randomize the initial configuration before sampling, of only the molecules
# we are interested in (Rpb4 and Rpb7)
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    excluded_rigid_bodies=fixed_rbs,
                                    max_translation=50,
                                    verbose=False,
                                    cutoff=5.0,
                                    niterations=100)

outputobjects = []  # reporter objects (for stat files)

# -----------------------------------
# Define Scoring Function Components
# -----------------------------------

# Here we are defining a number of restraints on our system.
#  For all of them we call add_to_model() so they are incorporated
#  into scoring. We also add them to the outputobjects list, so they
#  are reported in stat files

# Connectivity keeps things connected along the backbone (ignores if inside
# same rigid body)
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
#  To speed up this expensive restraint, we operate it at resolution 20
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         included_objects=root_hier,
                                         resolution=10)
ev.add_to_model()
outputobjects.append(ev)

cldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Protein 1")
cldbkc.set_residue1_key("Residue 1")
cldbkc.set_protein2_key("Protein 2")
cldbkc.set_residue2_key("Residue 2")
#cldbkc.set_id_score_key("Score")

cl2db = IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cl2db.create_set_from_file(datadirectory+'pdb2cl.csv')

cl = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,
                                   database=cl2db,
                                   length=21.0,
                                   slope=0.01,
                                   resolution=1.0,
                                   label="pdb2cl",
                                   linker=None,
                                   weight=1.)
cl.add_to_model()
outputobjects.append(cl)

# --------------------------
# Monte-Carlo Sampling
# --------------------------

# --------------------------
# Set MC Sampling Parameters
# --------------------------
num_frames = 20000
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

