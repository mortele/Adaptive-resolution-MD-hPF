import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Add the src/ directory to the python path so we can import the code 
# we need to use directly as 'from <file name> import <function/class>'
sys.path.append(os.path.join(os.path.dirname(__file__)))

from positions_and_velocities_and_forces import compileAndRun, getValuesLammps, readValues


if __name__ == '__main__':
    compileAndRun()

    md_energies = readValues("energy.dump", shape=(2,10))
    lammps_energies = readValues("lammps.energy", shape=(3,11))
    lammps_energies = lammps_energies[1:,1:]

    markersize = 3.0
    low_limit = 0.9e-16
    up_limit = 1e-0
    fontsize = 16

    plt.rc('text', usetex=True)
    plt.figure()
    difference = np.abs(md_energies[0,:] - lammps_energies[0,:])
    difference[difference<1e-16] = 1e-16
    plt.semilogy(difference, 'b-o', markersize=markersize)
    plt.ylabel(r"$|T_{md} - T_{lammps}|$", fontsize=fontsize)
    plt.subplots_adjust(left=0.3)
    plt.ylim((low_limit, up_limit))
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.savefig(os.path.join(os.path.dirname(__file__), '..',  'figures', 'kinetic-energy.png'), transparent=True, bbox_inches='tight')
    #plt.show()
    
    plt.figure()
    difference = np.abs(md_energies[1,:] - lammps_energies[1,:])
    difference[difference<1e-16] = 1e-16
    plt.semilogy(difference, 'b-o', markersize=markersize)
    plt.ylabel(r"$|V_{md} - V_{lammps}|$", fontsize=fontsize)
    plt.subplots_adjust(left=0.3)
    plt.ylim((low_limit, up_limit))
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.savefig(os.path.join(os.path.dirname(__file__),  '..', 'figures', 'potential-energy.png'), transparent=True, bbox_inches='tight')


    lammps_total = np.sum(lammps_energies, axis=0) 
    md_total     = np.sum(md_energies, axis=0)
    md_diff      = np.array([np.abs(md_total[i]-md_total[0]) for i in range(len(md_total))])
    md_diff[0]   = 1e-16
    lammps_diff   = np.array([np.abs(lammps_total[i]-lammps_total[0]) for i in range(len(lammps_total))])
    lammps_diff[0] = 1e-16

    print(lammps_total)

    plt.figure()
    difference = np.abs(md_total - lammps_total)
    difference[difference<1e-16] = 1e-16
    plt.semilogy(difference, 'b-o', markersize=markersize, label=r'$|E_{md}-E_{lammps}|$')
    plt.semilogy(md_diff,    'r-o', markersize=markersize, label=r'$|\Delta E_{md}|$')
    plt.semilogy(lammps_diff,    'k-o', markersize=markersize, label=r'$|\Delta E_{lammps}|$')
    plt.ylabel(r"$|\Delta E_{total}|$", fontsize=fontsize)
    plt.subplots_adjust(left=0.3)
    plt.ylim((low_limit, up_limit))
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.legend(fontsize=fontsize)
    plt.savefig(os.path.join(os.path.dirname(__file__),  '..', 'figures', 'total-energy.png'), transparent=True, bbox_inches='tight')
    plt.show()
