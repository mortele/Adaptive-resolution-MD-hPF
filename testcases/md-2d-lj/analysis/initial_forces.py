import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Compile test case and run LAMMPS to get output
os.chdir("../")
lmp = "/usr/local/bin/lmp_serial"
subprocess.run("make clean",         shell=True, stdout=subprocess.PIPE)
subprocess.run("make",               shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
subprocess.run(lmp + " < lammps.in", shell=True)
subprocess.run("./md-2d-lj.app",     shell=True, stdout=subprocess.PIPE)
os.chdir("analysis")

# Read positions/forces.dump and lammps.dump 
forces_md        = np.empty(shape=(3,100))
forces_lammps    = np.empty(shape=(3,100))
positions_md     = np.empty(shape=(3,100))
positions_lammps = np.empty(shape=(3,100))

with open(os.path.join(os.path.dirname(__file__), '..', 'forces.dump')) as inFile :
    content = inFile.readlines()
    for i in range(len(content)) :
        forces_md[:,i] = np.array([float(s) for s in content[i].split()])

with open(os.path.join(os.path.dirname(__file__), '..', 'positions.dump')) as inFile :
    content = inFile.readlines()
    for i in range(len(content)) :
        positions_md[:,i] = np.array([float(s) for s in content[i].split()])


with open(os.path.join(os.path.dirname(__file__), '..', 'lammps.dump')) as inFile :
    content = inFile.readlines()
    content = content[9:]
    content = content[100+9:]
    for i in range(100) :
        line                    = np.array([float(s) for s in content[i].split()])
        positions_lammps[:,i]   = line[3:6]
        forces_lammps[:,i]      = line[9:12]


plt.figure()
difference = np.abs(np.ravel(positions_md)-np.ravel(positions_lammps))
difference[difference<1e-16] = 1e-16
plt.semilogy(difference, 'r.')
plt.ylabel("|x(md) - x(lammps)|")
plt.subplots_adjust(left=0.2)

plt.figure()
difference = np.abs(np.ravel(forces_md)-np.ravel(forces_lammps))
difference[difference<1e-16] = 1e-16
plt.semilogy(difference, 'r.')
plt.ylabel("|F(md) - F(lammps)|")
plt.subplots_adjust(left=0.2)
plt.show()