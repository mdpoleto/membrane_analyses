import simtk.openmm as mm
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
import mdtraj as md

import os, argparse, sys
import numpy as np


ap = argparse.ArgumentParser(description=__doc__)

ap.add_argument('-crd', type=str, default=None, required=True,
                help='Input coordinate file (.crd)')
ap.add_argument('-psf', type=str, default=None, required=True,
                help='Topology file in XPLOR format (.psf)')
ap.add_argument('-dcd', type=str, default=None, required=True,
                help='Trajectory file in DCD format (.dcd)')
ap.add_argument('-toppar', type=str, default='toppar.str', required=True,
                help='Force field stream file (ex. "toppar.str").')
ap.add_argument('-boxsize', type=float, required=True, nargs=3,
                help='Box vector lengths (a,b,c) in nanometers')
ap.add_argument('-mol', type=str, default=None, required=True,
                help='Resname of your molecule of interest.')

ap.add_argument('-which_forces', type=str, default='all', required=True,
                help='All forces or just Nonbonded forces')

# Options
ap.add_argument('-temp', default=298, type=float,
                help='Target temperature, in Kelvin. Default = 298.')
ap.add_argument('-pressure', default=1.0, type=float,
                help='Target pressure, in bar. Default is 1.0.')
ap.add_argument('-dt', default=2, type=int,
                help='Integration step (in fs). Default = 1.')
ap.add_argument('-is_drude', default=False, type=str,
                help='If it is a Drude system or not')

cmd = ap.parse_args()

########################################################
def gen_box(psf, boxsize):

	boxvectors = boxsize

	boxlx = boxvectors[0]*nanometer
	boxly = boxvectors[1]*nanometer
	boxlz = boxvectors[2]*nanometer

	psf.setBox(boxlx, boxly, boxlz)
	return psf

def read_toppar(filename):
	extlist = ['rtf', 'prm', 'str']

	parFiles = ()
	for line in open(filename, 'r'):
		if '!' in line: line = line.split('!')[0]
		parfile = line.strip()
		if len(parfile) != 0:
			ext = parfile.lower().split('.')[-1]
			if not ext in extlist: continue
			parFiles += ( parfile, )

	params = CharmmParameterSet( *parFiles )
	return params, parFiles

def forcegroupify(system):
    forcegroups = {}
    for i, force in enumerate(system.getForces()):
            force.setForceGroup(i)
            #print(force.getName(), i)
            forcegroups[i] = force.getName()
    return forcegroups

def add_to_dict_vstack(dict, key, array):
    """ Insert an array into a dictionary accounting for existing keys """
    if key in dict.keys():
        prev_array = dict[key]
        new_array = np.vstack((prev_array, array))
        dict[key] = new_array
    else:
        dict[key] = array

    return dict

def t_or_f(arg):
    ua = str(arg).upper()
    if 'TRUE'.startswith(ua):
       return True
    elif 'FALSE'.startswith(ua):
       return False
    else:
       sys.exit("Flag -is_drude only accepts True or False.")

#########################################################

psffile = cmd.psf
crdfile = cmd.crd
dcdfile = cmd.dcd

temperature = cmd.temp*kelvin
pressure	= cmd.pressure*bar
dt          = cmd.dt
is_drude    = t_or_f(cmd.is_drude)
molecule    = cmd.mol


if cmd.which_forces.lower() == 'all' or cmd.which_forces.lower() == 'nonbonded':
    calc_mode = cmd.which_forces.lower()
else:
    sys.exit("Calculation mode not valid. Flag -which_forces must be either 'all' or 'nonbonded'.")
########################################################
print("\n> Setting the system:\n")
print(" -> Reading force field directory...")
charmm_params = read_toppar(cmd.toppar)[0]

# reading and parsing topology and coordinates
psf = CharmmPsfFile(psffile)
crd = CharmmCrdFile(crdfile)

print(" -> Setting box (using user's information)...")
psf = gen_box(psf, cmd.boxsize)

print(" -> Creating system and setting parameters...")
system = psf.createSystem(charmm_params, nonbondedMethod=PME, nonbondedCutoff=1.2*nanometers, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0005, constraints=HBonds)

# making sure each ForceGroup has its own index
forcegroups = forcegroupify(system)

if is_drude:
    integrator = DrudeLangevinIntegrator(temperature, 5/picosecond, 1*kelvin, 20/picosecond, 0.001*picoseconds)
    integrator.setMaxDrudeDistance(0.02) # Drude Hardwall
else:
    integrator = LangevinIntegrator(temperature, 1/picosecond, dt)

simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(crd.positions)

######################################
# Calculating Energies for entire trajectory
traj = md.load(cmd.dcd, top=cmd.psf)
top = traj.topology
ligand = top.select(str(molecule))

print("\n -> Atom IDs of selection:", ligand)


frame_dict = {}

print("\n>>> Reading and Calculating:")
for frame in range(traj.n_frames):
    print("> Frame " + str(frame))

    id_dict = {}

    # update simulation context with new positions
    simulation.context.setPositions(traj.xyz[frame])

    # update simulation context with new box vectors
    a, b, c = traj.unitcell_vectors[frame]
    simulation.context.setPeriodicBoxVectors(*(a,b,c))

    # get updated energies
    #state = simulation.context.getState(getEnergy=True)
    #potenergy = getPotentialEnergy(state)

    #forces_to_be_included = ['CustomNonbondedForce','NonbondedForce','DrudeForce']
    forces_to_be_included = ["HarmonicBondForce",     # Bonds
                             "HarmonicAngleForce",    # Angles
                             #"HarmonicBondForce",    # Urey-Bradley, it is coverd by Bonds
                             "PeriodicTorsionForce",  # Dihedrals
                             "CustomTorsionForce",    # Impropers
                             "CMAPTorsionForce",      # CMAP
                             "NonbondedForce",        # Elec + VDW + some Drude forces
                             "DrudeForce",            # The remaining Drude forces
                             "CustomNonbondedForce"]  # NBFixes of any kinds"""


    if calc_mode == 'all':
        forces = simulation.context.getState(getForces=True).getForces(asNumpy=True)
        forces = forces.value_in_unit(kilocalorie/(angstrom*mole)) # remove unit

        for id in ligand:
            id_dict = add_to_dict_vstack(id_dict, id, forces[id])


    else: # if calc_mode == nonbonded
        tmp_force_contributions = {}
        for forceid, forcename in forcegroups.items():
            if forcename in forces_to_be_included:
                forces = simulation.context.getState(getForces=True, groups={forceid}).getForces(asNumpy=True)
                forces = forces.value_in_unit(kilocalorie/(angstrom*mole)) # remove unit
                #print(forceid, forcename)

                for id in ligand:
                    #print(top.atom(id).name, forces[id])
                    tmp_force_contributions = add_to_dict_vstack(tmp_force_contributions, id, forces[id])


        # Sum forces contributions
        for id, force_array in tmp_force_contributions.items():
            total_force = np.sum(force_array, axis=0)
            id_dict = add_to_dict_vstack(id_dict, id, total_force)

    frame_dict[frame] = id_dict

#print(frame_dict)
########################################################################

o = open("Forces_overview.dat", "w")

if is_drude:
    drude_particles = np.array([], dtype=int)
    for id in ligand:
        if top.atom(id).name[0] == "D":
            drude_particles = np.append(drude_particles, id)

    ligand = np.setdiff1d(ligand, drude_particles)

print("\n###############################")
print("NAME,                  AVG,                        STD")
for id in ligand:
    tmp_avg = np.empty([0,4])

    atomname = top.atom(id).name
    resname = top.atom(id).residue.name
    globalresid = top.atom(id).residue.index + 1

    identifier = resname + "_" + str(globalresid) + "_" + atomname
    ot = open(identifier + ".dat", "w")
    ot.write("#Frame            FR             FX             FY             FZ\n")

    for frame, id_dict in frame_dict.items():

        if not is_drude:
            fx,fy,fz = id_dict[id]

        else: # if it is a Drude system
            drude_id = id + 1 # hard-coded feature of CHARMM: Drude particles
                              # are always built after their parent atom, so
                              # Drude_ID = Parent_ID + 1 always
            if drude_id in drude_particles:
                fxp,fyp,fzp = id_dict[id]
                fxd,fyd,fzd = id_dict[id+1]
                fx = fxp + fxd
                fy = fyp + fyd
                fz = fzp + fzd
            else:
                fx,fy,fz = id_dict[id]

        fr = np.sqrt(fx**2 + fy**2 + fz**2)
        fr = round(fr,4)
        string = "{a:<12}{b: 15.4f}{c: 15.4f}{d: 15.4f}{e: 15.4f}\n".format(a=frame, b=fr, c=fx, d=fy,e=fz)
        ot.write(string)

        tmp_avg = np.vstack((tmp_avg, np.array([fr, fx, fy, fz])))

    ot.close()

    avg = np.mean(tmp_avg, axis=0)
    std = np.std(tmp_avg,  axis=0)

    string = identifier + "  ,  " + str(avg) + "  ,  " + str(std)
    print(string)
    o.write(string + "\n")
