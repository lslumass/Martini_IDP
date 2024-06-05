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
chk_file = 'md.chk'

# simulation setups
nvt_step = 50000
npt_step = 50000
prod_step = 250000000
log_freq = 1000
xtc_freq = 5000                     # output frequency
chk_freq = 5000000                  # checkpoint file frequency

dt = 20*femtosecond                 # time step
epsilon_r = 15                      # relative electric constant
T = 303*kelvin                      # system temperature
friction = 0.1/picosecond           # friction constant for LangevinIntegrator

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
integrator = LangevinIntegrator(T, friction, dt)
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
sim.context.setVelocitiesToTemperature(T)
print('Running NVT equilibration...')
sim.step(nvt_step)

################################################################################
### NPT equilibration ###
### select bilayer or non-bilayer system ###
# for bilayer system
#barostat = MonteCarloMembraneBarostat(1 * bar, 0 * bar * nanometer, T,
#		MonteCarloMembraneBarostat.XYIsotropic,
#		MonteCarloMembraneBarostat.ZFree, 10
#		)
# for non-bilayer system
barostat = MonteCarloBarostat(1 * bar, T)

system.addForce(barostat)
sim.context.reinitialize(True)
print('Running NPT equilibration...')
sim.step(npt_step)
sim.saveState('equi.state')
sim.saveCheckpoint('equi.chk')

################################################################################
### Production run ###
sim.reporters.append(StateDataReporter(log_file, log_freq, progress=True, step=True, potentialEnergy=True, totalEnergy=True, temperature=True, speed=True, totalSteps=prod_step))
sim.reporters.append(XTCReporter(xtc_file, xtc_freq))
sim.reporters.append(CheckpointReporter(chk_file, chk_freq))

print("Running simulation...")
sim.step(prod_step)

sim.saveCheckpoint('final.chk')
print('Finish!')