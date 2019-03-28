import logging
import os
from pdbfixer import PDBFixer
import simtk.openmm as mm
from simtk.openmm import unit
from simtk.openmm import version
from simtk.openmm.app import Topology, PDBFile, Modeller, ForceField, PDBxFile, PME, Simulation, StateDataReporter
from openmmtools import states, mcmc
import protein_features as pf
import matplotlib.pyplot as plot
import numpy as np
import tempfile

## Setup general logging (guarantee output/error message in case of interruption)
logger = logging.getLogger(__name__)
logging.root.setLevel(logging.DEBUG)
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("urllib3").setLevel(logging.WARNING)

'''
## clean up the input pdb file using pdbfixer and load using Modeller

# fix using pdbfixer: remove the ligand but keep the crystal waters
fixer = PDBFixer(filename='3UG9_clean.pdb')
fixer.findMissingResidues()

# modify missingResidues so the extra residues on the end are ignored
#fixer.missingResidues = {(0,47): fixer.missingResidues[(0,47)]}
fixer.missingResidues = {}

# remove ligand but keep crystal waters
fixer.removeHeterogens(True)
print("Done removing heterogens.")

# find missing atoms/terminals
fixer.findMissingAtoms()
if fixer.missingAtoms or fixer.missingTerminals:
    fixer.addMissingAtoms()
    print("Done adding atoms/terminals.")
else:
    print("No atom/terminal needs to be added.")

# add hydrogens
fixer.addMissingHydrogens(7.0)
print("Done adding hydrogens.")
# output fixed pdb
PDBFile.writeFile(fixer.topology, fixer.positions, open('5UG9_fixed.pdb', 'w'), keepIds=True)
print("Done outputing the fixed pdb file.")

# load pdb to Modeller
pdb = PDBFile('5UG9_fixed.pdb')
molecule = Modeller(pdb.topology,pdb.positions)
print("Done loading pdb to Modeller.")
# load force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
print("Done loading force field.")
print("OpenMM version:", version.version)
# prepare system
molecule.addSolvent(forcefield, padding=12*unit.angstrom, model='tip3p', positiveIon='Na+', negativeIon='Cl-', ionicStrength=0*unit.molar)
print("Done adding solvent.")
PDBxFile.writeFile(molecule.topology,molecule.positions,open("5UG9_fixed.pdbx", 'w'))
PDBFile.writeFile(molecule.topology,molecule.positions,open("5UG9_fixed_solvated.pdb", 'w'))
print("Done outputing pdbx and solvated pdb.")
system = forcefield.createSystem(molecule.topology, nonbondedMethod=PME, rigidWater=True, nonbondedCutoff=1*unit.nanometer)

# specify the rest of the context for minimization
integrator = mm.VerletIntegrator(0.5*unit.femtoseconds)
print("Done specifying integrator.")
platform = mm.Platform.getPlatformByName('CUDA')
print("Done specifying platform.")
platform.setPropertyDefaultValue('Precision', 'single')
print("Done setting the precision to single.")
minimize = Simulation(molecule.topology, system, integrator, platform)
print("Done specifying simulation.")
minimize.context.setPositions(molecule.positions)
print("Done recording a context for positions.")
minimize.context.setVelocitiesToTemperature(310.15*unit.kelvin)
print("Done assigning velocities.")

# start minimization
tolerance = 0.1*unit.kilojoules_per_mole/unit.angstroms
print("Done setting tolerance.")
minimize.minimizeEnergy(tolerance=tolerance,maxIterations=1000)
print("Done setting energy minimization.")
minimize.reporters.append(StateDataReporter('relax-hydrogens.log', 1000, step=True, temperature=True, potentialEnergy=True, totalEnergy=True, speed=True))
minimize.step(100000)
print("Done 100000 steps of minimization.")
print("Potential energy after minimization:")
print(minimize.context.getState(getEnergy=True).getPotentialEnergy())
positions = minimize.context.getState(getPositions=True).getPositions()
print("Done updating positions.")
#velocities = minimize.context.getState(getVelocities=True).getVelocities()
#print("Done updating velocities.")
minimize.saveCheckpoint('state.chk')
print("Done saving checkpoints.")
# update the current context with changes in system
#minimize.context.reinitialize(preserveState=True)
# output the minimized protein as a shortcut
PDBFile.writeFile(molecule.topology,positions,open("5UG9_minimized.pdb", 'w'))
print("Done outputing minimized pdb.")
'''
# directly load the minimized protein
pdb = PDBFile('5UG9_minimized.pdb')
molecule = Modeller(pdb.topology,pdb.positions)
# load force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
print("Done loading force field.")
print("OpenMM version:", version.version)
system = forcefield.createSystem(molecule.topology, nonbondedMethod=PME, rigidWater=True, nonbondedCutoff=1*unit.nanometer)

# specify the collective variables
temperature = 310.15 * unit.kelvin
pressure = 1.0 * unit.atmospheres

expt = 'dih1'
min_v = -178.7/180*np.pi #this is the estimated min value in radians
max_v = -56.6/180*np.pi # also estimated
sigma = 14.79/180*np.pi
nstates = 50
iteration = 1000
# add the custom cv force
# Specify the set of key atoms and calculate key dihedrals and distances
(dih, dis) = pf.main('5UG9','A')
# dihedrals
kT = unit.MOLAR_GAS_CONSTANT_R * temperature
kT_in_md_units = kT.value_in_unit_system(unit.md_unit_system)

# Create Dunbrack biasing potential
# Parameters for each torsion are
# 'dih0_phi0' - torsion center of state 0 (BLAminus) for dih0
# 'dih0_dphi' - stddev of state 0 (BLAminus) for dih0
# ...

initial_state = 0

# TODO:
# dunbrack_phi0[state_index,dihedral_index] is the reference dihedral value for state 'state_index' and dihedral 'dihedral_index'
# dunbrack_dphi[state_index,dihedral_index] is the dihedral stddev for state 'state_index' and dihedral 'dihedral_index'
# 'ndihedrals' is the number of dihedrals we want to restrain
# dihedral_atoms[dihedral_index] is 4-tuple or list of atom indices involved in dihedral 'dihedral_index'


for dihedral_index in range(ndihedrals):
    energy_expression = '-lambda*(K/2)*cos(theta-dih{dihedral_index}_phi0); K={kT_in_md_units}/dih{dihedral_index}_dphi**2'.format(vars())
    force = mm.CustomTorsionForce(energy_expression)
    torsion_atom_indices = [ int(index) for inex in dihedral_atoms[dihedral_index] ]
    force.addTorsion(*torsion_atom_indices, [])
    force.addGlobalParameter('dih{dihedral_index}_phi0'.format(vars()), dunbrack_phi0[initial_state, dihedral_index])
    force.addGlobalParameter('dih{dihedral_index}_dphi'.format(vars()), dunbrack_dphi[initial_state, dihedral_index])
    force.addGlobalParameter('lambda', 0.0) # controls overall restraint strength
    system.addForce(force)
    print("Done defining the CV force.")

# create thermodynamic states
class MyComposableState(states.GlobalParameterState):
    rc_lambda = states.GlobalParameterState.GlobalParameter('lambda', standard_value=0.0)
    # TODO: Have Andrea check that this code below is valid for defining all the global parameters we want
    rc_dunbrack_phi0 = list()
    rc_dunbrack_dphi = list()
    for dihedral_index in range(ndihedrals):
        rc_dunbrack_phi0.append(states.GlobalParameterState.GlobalParameter('dih{dihedral_index}_phi0'.format(vars()), standard_value=0.0))
        rc_dunbrack_dphi.append(states.GlobalParameterState.GlobalParameter('dih{dihedral_index}_dphi'.format(vars()), standard_value=1.0))

protocol = dict()
# Add restraints to each of the Dunbrack states
for dihedral_index in range(ndihedrals):
    protocol['dih{dihedral_index}_phi0'.format(vars())] = list(dunbrack_phi0[:,dihedral_index])
    protocol['dih{dihedral_index}_dphi'.format(vars())] = list(dunbrack_dphi[:,dihedral_index])
protocol['lambda'] = np.ones([dunbrack_phi0.shape[0]])
# Add one more state where we turn off restraints
for dihedral_index in range(ndihedrals):
    protocol['dih{dihedral_index}_phi0'.format(vars())].append(0.0)
    protocol['dih{dihedral_index}_dphi'.format(vars())].append(1.0)
    protocol['lambda'].append(0.0) # turn off restraints

constants = {
    'temperature' : temperature,
    'pressure' : pressure,
}

composable_state = MyComposableState.from_system(system)
thermo_states = states.create_thermodynamic_state_protocol(system=system, protocol=protocol, constants=constants, composable_states=[composable_state])
# assign sampler_state
sampler_state = states.SamplerState(positions=pdb.positions, box_vectors=system.getDefaultPeriodicBoxVectors())

# Set up the context for mtd simulation
# at this step the CV and the system are separately passed to Metadynamics
from yank.multistate import SAMSSampler, MultiStateReporter
# TODO: You can use heavy hydrogens and 4 fs timesteps
move = mcmc.LangevinDynamicsMove(timestep=2.0*unit.femtoseconds, collision_rate=1.0/unit.picosecond, n_steps=1000, reassign_velocities=False)
print("Done specifying integrator for simulation.")
simulation = SAMSSampler(mcmc_moves=move, number_of_iterations=iteration, online_analysis_interval=None, gamma0=1.0, flatness_threshold=0.2)
storage_path = '/home/guoj1/projects/CV-selection/sams_simulation/5UG9_{}_{}/traj.nc'.format(expt,iteration)
reporter = MultiStateReporter(storage_path, checkpoint_interval=1)
# We should also add unsampled_states with the fully-interacting system
simulation.create(thermodynamic_states=thermo_states, sampler_states=[sampler_state], storage=reporter)
print("Done specifying simulation.")
print("Potential energy right before simulation:")
print(simulation.sampler_states[0].potential_energy)

# Run small-scale simulation (1000 iterations) and plot the free energy landscape
simulation.run()
print("Done with {} iterations.".format(iteration))
