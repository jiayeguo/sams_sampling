import logging
import os
from pdbfixer import PDBFixer
import simtk.openmm as mm
from simtk.openmm import unit 
from simtk.openmm import version
from simtk.openmm.app import Topology, PDBFile, Modeller, ForceField, PDBxFile, PME, Simulation, StateDataReporter, DCDReporter
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
fixer = PDBFixer(filename='5UG9_A.pdb')
fixer.findMissingResidues()

# modify missingResidues so the extra residues on the end are ignored
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
'''

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

temperature = 310.15 * unit.kelvin
# add the custom cv force
# Specify the set of key atoms and calculate key dihedrals and distances
(dih, dis) = pf.main('5UG9','A')
# dihedrals
kT = unit.MOLAR_GAS_CONSTANT_R * temperature 

K0 = kT/(-1.807 - (-1.949))**2
dih_0 = mm.CustomTorsionForce('(K0/2)*dphi ; dphi=cos(theta-tdih0)')
dih_0.addTorsion(int(dih[0][0]), int(dih[0][1]), int(dih[0][2]), int(dih[0][3]))
dih_0.addGlobalParameter('K0',K0)
dih_0.addGlobalParameter('tdih0',-1.949)

K1 = kT/(-2.246 - (-1.292))**2
dih_1 = mm.CustomTorsionForce('(K1/2)*dphi ; dphi=cos(theta-tdih1)')
dih_1.addTorsion(int(dih[1][0]), int(dih[1][1]), int(dih[1][2]), int(dih[1][3]))
dih_1.addGlobalParameter('K1',K1)
dih_1.addGlobalParameter('tdih1',-1.292)

K2 = kT/(2.404 - 2.403)**2
dih_2 = mm.CustomTorsionForce('(K2/2)*dphi ; dphi=cos(theta-tdih2)')
dih_2.addTorsion(int(dih[2][0]), int(dih[2][1]), int(dih[2][2]), int(dih[2][3]))
dih_2.addGlobalParameter('K2',K2)
dih_2.addGlobalParameter('tdih2',2.403)

K3 = kT/(-2.306 - (-1.962))**2
dih_3 = mm.CustomTorsionForce('(K3/2)*dphi ; dphi=cos(theta-tdih3)')
dih_3.addTorsion(int(dih[3][0]), int(dih[3][1]), int(dih[3][2]), int(dih[3][3]))
dih_3.addGlobalParameter('K3',K3)
dih_3.addGlobalParameter('tdih3',-1.962)
'''
dih_4 = mm.CustomTorsionForce("theta")
dih_4.addTorsion(int(dih[4][0]), int(dih[4][1]), int(dih[4][2]), int(dih[4][3]))

dih_5 = mm.CustomTorsionForce("theta")
dih_5.addTorsion(int(dih[5][0]), int(dih[5][1]), int(dih[5][2]), int(dih[5][3]))

dih_6 = mm.CustomTorsionForce("theta")
dih_6.addTorsion(int(dih[6][0]), int(dih[6][1]), int(dih[6][2]), int(dih[6][3]))

dih_7 = mm.CustomTorsionForce("theta")
dih_7.addTorsion(int(dih[7][0]), int(dih[7][1]), int(dih[7][2]), int(dih[7][3]))

# distances
dis_0 = mm.CustomBondForce("r")
dis_0.addBond(int(dis[0][0]), int(dis[0][1]))
dis_1 = mm.CustomBondForce("r")
dis_1.addBond(int(dis[1][0]), int(dis[1][1]))
dis_2 = mm.CustomBondForce("r")
dis_2.addBond(int(dis[2][0]), int(dis[2][1]))
dis_3 = mm.CustomBondForce("r")
dis_3.addBond(int(dis[3][0]), int(dis[3][1]))
dis_4 = mm.CustomBondForce("r")
dis_4.addBond(int(dis[4][0]), int(dis[4][1]))
print("Done populating dihedrals and distances.")
'''
# Specify a unique CustomCVForce
cv_force = mm.CustomCVForce('dih_0 + dih_1 + dih_2 + dih_3')
cv_force.addCollectiveVariable('dih_0', dih_0)
cv_force.addCollectiveVariable('dih_1', dih_1)
cv_force.addCollectiveVariable('dih_2', dih_2)
cv_force.addCollectiveVariable('dih_3', dih_3)
#cv_force.addCollectiveVariable('dih_4', dih_4)
#cv_force.addCollectiveVariable('dih_5', dih_5)
#cv_force.addCollectiveVariable('dih_6', dih_6)
#cv_force.addCollectiveVariable('dih_7', dih_7)
#cv_force.addCollectiveVariable('dis_0', dis_0)
#cv_force.addCollectiveVariable('dis_1', dis_1)
#cv_force.addCollectiveVariable('dis_2', dis_2)
#cv_force.addCollectiveVariable('dis_3', dis_3)
#cv_force.addCollectiveVariable('dis_4', dis_4)
print("Done defining forces.")

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
minimize.step(10000)
print("Done 10000 steps of simulation.")
#positions = minimize.context.getState(getPositions=True).getPositions()
#print("Done updating positions.")
#velocities = minimize.context.getState(getVelocities=True).getVelocities()
#print("Done updating velocities.") 
minimize.saveCheckpoint('state.chk')
print("Done saving checkpoints.")
system.addForce(cv_force)
print("Done adding CV force.")
# update the current context with changes in system
#minimize.context.reinitialize()
system = minimize.context.getSystem()
print("whether added:",system.getNumForces(),system.getForce(0),system.getForce(1),system.getForce(2),system.getForce(3),system.getForce(4),system.getForce(5))
#DEBUG
print("potential energy",minimize.context.getState(getEnergy=True, groups=32).getPotentialEnergy())

# create thermodynamic states
thermo_state = states.ThermodynamicState(system=system, temperature=temperature)
thermodynamic_states = []
thermodynamic_states.append(thermo_state)

# assign sampler_state
box_vectors = molecule.topology.getPeriodicBoxVectors()
#sampler_state = states.SamplerState(positions= positions, velocities = velocities, box_vectors = box_vectors)
sampler_state = states.SamplerState.from_context(minimize.context, ignore_collective_variables=False)
# Set up the context for mtd simulation
# at this step the CV and the system are separately passed to Metadynamics
from yank.multistate import SAMSSampler, MultiStateReporter
move = mcmc.LangevinDynamicsMove(timestep=0.002*unit.picoseconds, collision_rate=1.0/unit.picosecond, n_steps=1000, reassign_velocities=False)
print("Done specifying integrator.")
simulation = SAMSSampler(mcmc_moves=move, number_of_iterations=100000, online_analysis_interval=None, gamma0=1.0, flatness_threshold=0.2)
storage_path = tempfile.NamedTemporaryFile(delete=False).name + '.nc'
reporter = MultiStateReporter(storage_path, checkpoint_interval=1)
simulation.create(thermodynamic_states=thermodynamic_states, sampler_states=[sampler_state], storage=reporter)
print("Done specifying simulation.")

# Run small-scale simulation (10ns, 10^5 steps) and plot the free energy landscape
simulation.run()
print("Done with 10^5 steps of production run.")

