import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Compile test case and run LAMMPS to get output
os.chdir("../")
lmp = "/usr/local/bin/lmp_serial"
subprocess.run("make clean",         shell=True, stdout=subprocess.PIPE)
subprocess.run("make -j8",           shell=True, stdout=subprocess.PIPE)#, stderr=subprocess.PIPE)
subprocess.run(lmp + " < lammps.in", shell=True)
subprocess.run("./md-2d-lj.app",     shell=True, stdout=subprocess.PIPE)
os.chdir("analysis")

def getValuesLammps(step=0) :
    positions   = np.empty(shape=(3,100))
    velocities  = np.empty(shape=(3,100))
    forces      = np.empty(shape=(3,100))
    with open(os.path.join(os.path.dirname(__file__), '..', 'lammps.dump')) as inFile :
        content = inFile.readlines()
        content = content[9:]
        content = content[(100+9)*(step+1):]
        for i in range(100) :
            line            = np.array([float(s) for s in content[i].split()])
            positions[:,i]  = line[3:6]
            velocities[:,i] = line[6:9]
            forces[:,i]     = line[9:12]
    return positions, velocities, forces



# Read positions/forces.dump and lammps.dump 
forces_md        = np.empty(shape=(3,100))
positions_md     = np.empty(shape=(3,100))
positions2_md    = np.empty(shape=(3,100))
positions_lammps, _, forces_lammps = getValuesLammps(0)

with open(os.path.join(os.path.dirname(__file__), '..', 'forces.dump')) as inFile :
    content = inFile.readlines()
    for i in range(len(content)) :
        forces_md[:,i] = np.array([float(s) for s in content[i].split()])

with open(os.path.join(os.path.dirname(__file__), '..', 'positions.dump')) as inFile :
    content = inFile.readlines()
    for i in range(len(content)) :
        positions_md[:,i] = np.array([float(s) for s in content[i].split()])

with open(os.path.join(os.path.dirname(__file__), '..', 'positions2.dump')) as inFile :
    content = inFile.readlines()
    for i in range(len(content)) :
        positions2_md[:,i] = np.array([float(s) for s in content[i].split()])


plt.figure()
difference = np.abs(np.ravel(forces_md)-np.ravel(forces_lammps))
difference[difference<1e-16] = 1e-16
plt.semilogy(difference, 'r.')
plt.ylabel("|F(md) - F(lammps)|")
plt.subplots_adjust(left=0.3)

plt.figure()
difference = np.abs(np.ravel(positions_md)-np.ravel(positions_lammps))
difference[difference<1e-16] = 1e-16
plt.semilogy(difference, 'r.')
plt.ylabel("|x(md) - x(lammps)|")
plt.subplots_adjust(left=0.3)

plt.figure()
positions_lammps, _, _ = getValuesLammps(1)
difference = np.abs(np.ravel(positions2_md)-np.ravel(positions_lammps))
difference[difference<1e-16] = 1e-16
plt.semilogy(difference, 'r.')
plt.ylabel("|x2(md) - x2(lammps)|")
plt.subplots_adjust(left=0.3)
plt.show()