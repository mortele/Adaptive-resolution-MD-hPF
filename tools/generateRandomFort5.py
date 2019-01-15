import numpy as np

# Parameters
n_molecules          = 10
n_atoms_per_molecule = 1
atom_label           = 'A'
atom_id              = 1
n_bonds              = 0
max_n_bonds          = 6
box_size = np.array([10.0, 10.0, 10.0])

print("# Simulation box: x, y, z, binc")
for i in range(3) :
    print(box_size[i],end=" ")
print("0.0")

print("# Total number of molecules: N")
print(n_molecules)

atom_counter = 1
for molecule in range(1, n_molecules+1) :
    print("# Molecule number #", molecule)
    print(n_atoms_per_molecule)

    print(atom_counter, end=" ")
    atom_counter += 1

    print(atom_label, end=" ")
    print(atom_id,    end=" ")
    print(n_bonds,    end=" ")

    for i in range(3) :
        coordinate = np.random.uniform(low=0.0, high=box_size[i])
        print(coordinate, end=" ")
    for i in range(max_n_bonds) :
        print("0 ", end="")
    print("")



