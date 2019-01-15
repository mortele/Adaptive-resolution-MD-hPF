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

coordinates = np.random.uniform(low=0.0, high=10.0, size=(3, n_molecules))
atom_counter = 1
for molecule in range(1, n_molecules+1) :
    print("# Molecule number #", molecule)
    print(n_atoms_per_molecule)

    print(atom_counter, end=" ")
    print(atom_label, end=" ")
    print(atom_id,    end=" ")
    print(n_bonds,    end=" ")

    # x,y,z
    for i in range(3) :
        print(coordinates[i,atom_counter-1], end=" ")

    # vx, vy, vz
    for i in range(3) :
        print("0.0 ", end=" ")

    # bonds, all zeros--has to be exactly max_n_bonds for each atom
    for i in range(max_n_bonds) :
        print("0 ", end="")
    print("")
    atom_counter += 1


print("\n\n\n")

for i in range(n_molecules) :
    print("positions(:,",i+1,") = [", end="")
    for j in range(2) :
        print(coordinates[j,i], end="_real64,")
    print(coordinates[-1,i], end="_real64 ]\n")

positions(:, 1 )  = [ 9.861513516680217_real64,  5.854609220329756_real64,  4.090152684461406_real64   ]
positions(:, 2 )  = [ 2.9071826059146924_real64, 8.743543693640504_real64,  9.004829883778097_real64   ]
positions(:, 3 )  = [ 2.461251549282577_real64,  6.473427427545419_real64,  0.40419124737701484_real64 ]
positions(:, 4 )  = [ 2.873699847732466_real64,  4.562054499100865_real64,  0.8750826555311375_real64  ]
positions(:, 5 )  = [ 8.61517700020657_real64,   7.5947341584529235_real64, 7.817588697757908_real64   ]
positions(:, 6 )  = [ 8.329323455469527_real64,  6.48491076670158_real64,   1.0595951799017522_real64  ]
positions(:, 7 )  = [ 7.539162628881729_real64,  6.8390933240322624_real64, 7.783573997831711_real64   ]
positions(:, 8 )  = [ 5.502628575319229_real64,  5.086691204907025_real64,  9.508740311058814_real64   ]
positions(:, 9 )  = [ 4.68875984367639_real64,   0.7111235962240725_real64, 7.499310540550107_real64   ]
positions(:, 10 ) = [ 5.628789714994156_real64,  0.5249494827666557_real64, 6.096779978340795_real64   ]




