from openmm.unit import *
from openmm import *
from openmm.app import *
import martini_openmm as martini
from mdtraj.reporters import XTCReporter
from sys import stdout


# file names
input_gro = 'system.gro'
top_file = 'topol.top'
log_file = 'md.log'
xtc_file = 'md.xtc'

nvt_step = 50000
npt_step = 50000
prod_step = 250000000
epsilon_r = 15
platform = Platform.getPlatformByName("CUDA")
properties = {'Precision': 'mixed'}
conf = GromacsGroFile(input_gro)
box_vectors = conf.getPeriodicBoxVectors()

# get any defines
defines = {}
try:
    with open("defines.txt") as def_file:
        for line in def_file:
            line = line.strip()
            defines[line] = True
except FileNotFoundError:
    pass

top = martini.MartiniTopFile(top_file, periodicBoxVectors=box_vectors, defines=defines, epsilon_r=epsilon_r)
system = top.create_system(nonbonded_cutoff=1.1*nanometer)
integrator = LangevinIntegrator(303*kelvin, 0.1/picosecond, 20*femtosecond)
integrator.setRandomNumberSeed(0)
sim = Simulation(top.topology, system, integrator, platform, properties)
sim.context.setPositions(conf.getPositions())

################################################################################
### Minimization ###
#sim.reporters.append(PDBReporter('em.pdb', 1000))
sim.reporters.append(StateDataReporter(stdout, 5000, step=True, potentialEnergy=True, temperature=True, volume=True))
print("Minimizing energy...")
sim.minimizeEnergy(maxIterations=5000, tolerance=1.0)
energies = sim.context.getState(getEnergy=True).getPotentialEnergy()
print("System minimized at", energies, "\n")

################################################################################
### NVT equilibration ###
sim.context.setVelocitiesToTemperature(303*kelvin)
print('Running NVT equilibration...')
sim.step(nvt_step)

################################################################################
### NPT equilibration ###
barostat = openmm.MonteCarloMembraneBarostat(1 * bar, 0 * bar * nanometer, 303 * kelvin,
		openmm.MonteCarloMembraneBarostat.XYIsotropic,
		openmm.MonteCarloMembraneBarostat.ZFree, 10
		)
system.addForce(barostat)
sim.context.reinitialize(True)
print('Running NPT equilibration...')
sim.step(npt_step)
sim.saveState('equi.state')
sim.saveCheckpoint('equi.chk')

################################################################################
### Production run ###
sim.reporters.append(StateDataReporter(log_file, 1000, step=True, potentialEnergy=True, totalEnergy=True, density=True, temperature=True, volume=True))
xtc_reporter = XTCReporter(xtc_file, 100000)
sim.reporters.append(xtc_reporter)
print("Running simulation...")
sim.step(prod_step)
