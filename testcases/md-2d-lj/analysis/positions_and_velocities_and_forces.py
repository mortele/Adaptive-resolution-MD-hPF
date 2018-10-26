import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Compile test case and run LAMMPS to get output
def compileAndRun() :
    os.chdir("../")
    lmp = "/usr/local/bin/lmp_serial"
    subprocess.run("make clean",         shell=True, stdout=subprocess.PIPE)
    subprocess.run("make -j8",           shell=True, stdout=subprocess.PIPE)#, stderr=subprocess.PIPE)
    subprocess.run(lmp + " < lammps.in", shell=True)
    subprocess.run("./md-2d-lj.app",     shell=True)
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

def readValues(fileName, shape=(3,100)) :
    values = np.empty(shape=shape)
    with open(os.path.join(os.path.dirname(__file__), '..', fileName)) as inFile :
        content = inFile.readlines()
        j = 0
        for i in range(len(content)) :
            if content[i][0] != '#' :
                if content[i].strip() :
                    values[:,j] = np.array([float(s) for s in content[i].split()])
                    j += 1
    return values


if __name__ == '__main__':
    compileAndRun()

    # Read positions/forces.dump and lammps.dump 
    positions_md     = readValues('positions.dump')
    velocities_md    = readValues('velocities.dump')
    forces_md        = readValues('forces.dump')

    positions2_md    = readValues('positions2.dump')
    velocities2_md   = readValues('velocities2.dump')
    forces2_md       = readValues('forces2.dump')

    positions_lammps,  velocities_lammps,  forces_lammps  = getValuesLammps(1)
    positions2_lammps, velocities2_lammps, forces2_lammps = getValuesLammps(2)

    print(positions_lammps[:,0])
    markersize = 3.0
    low_limit = 0.9e-16
    up_limit = 1e-0
    fontsize = 16

    plt.rc('text', usetex=True)

    # POSITION STEP 1
    plt.figure()
    difference = np.abs(np.ravel(positions_md)-np.ravel(positions_lammps))
    difference[difference<1e-16] = 1e-16
    plt.semilogy(difference, 'r.', markersize=markersize)
    plt.ylabel(r"$|r_{1,md} - r_{1,lammps}|$", fontsize=fontsize)
    plt.subplots_adjust(left=0.3)
    plt.ylim((low_limit, up_limit))
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.savefig(os.path.join(os.path.dirname(__file__), 'figures', 'positions-firststep.png'), transparent=True, bbox_inches='tight')

    # VELOCITY STEP 1
    plt.figure()
    difference = np.abs(np.ravel(velocities_md)-np.ravel(velocities_lammps))
    difference[difference<1e-16] = 1e-16
    plt.semilogy(difference, 'r.', markersize=markersize)
    plt.ylabel(r"$|v_{1,md} - v_{1,lammps}|$", fontsize=fontsize)
    plt.subplots_adjust(left=0.3)
    plt.ylim((low_limit, up_limit))
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.savefig(os.path.join(os.path.dirname(__file__), 'figures', 'velocities-firststep.png'), transparent=True, bbox_inches='tight')
    
    # FORCES STEP 1
    plt.figure()
    difference = np.abs(np.ravel(forces_md)-np.ravel(forces_lammps))
    difference[difference<1e-16] = 1e-16
    plt.semilogy(difference, 'r.', markersize=markersize)
    plt.ylabel(r"$|F_{1,md} - F_{1,lammps}|$", fontsize=fontsize)
    plt.ylim((low_limit, up_limit))
    plt.subplots_adjust(left=0.3)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.savefig(os.path.join(os.path.dirname(__file__), 'figures', 'forces-firststep.png'), transparent=True, bbox_inches='tight')

    # POSITION STEP 2
    plt.figure()
    difference = np.abs(np.ravel(positions2_md)-np.ravel(positions2_lammps))
    difference[difference<1e-16] = 1e-16
    plt.semilogy(difference, 'r.', markersize=markersize)
    plt.ylabel(r"$|r_{2,md} - r_{2,lammps}|$", fontsize=fontsize)
    plt.subplots_adjust(left=0.3)
    plt.ylim((low_limit, up_limit))
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.savefig(os.path.join(os.path.dirname(__file__), 'figures', 'positions-secondstep.png'), transparent=True, bbox_inches='tight')

    # VELOCITY STEP 2
    plt.figure()
    difference = np.abs(np.ravel(velocities2_md)-np.ravel(velocities2_lammps))
    difference[difference<1e-16] = 1e-16
    plt.semilogy(difference, 'r.', markersize=markersize)
    plt.ylabel(r"$|v_{2,md} - v_{2,lammps}|$", fontsize=fontsize)
    plt.subplots_adjust(left=0.3)
    plt.ylim((low_limit, up_limit))
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.savefig(os.path.join(os.path.dirname(__file__), 'figures', 'velocities-secondstep.png'), transparent=True, bbox_inches='tight')

    # FORCES STEP 2
    plt.figure()
    difference = np.abs(np.ravel(forces2_md)-np.ravel(forces2_lammps))
    difference[difference<1e-16] = 1e-16
    plt.semilogy(difference, 'r.', markersize=markersize)
    plt.ylabel(r"$|F_{2,md} - F_{2,lammps}|$", fontsize=fontsize)
    plt.subplots_adjust(left=0.3)
    plt.ylim((low_limit, up_limit))
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.savefig(os.path.join(os.path.dirname(__file__), 'figures', 'forces-secondstep.png'), transparent=True, bbox_inches='tight')
    
    plt.show()

