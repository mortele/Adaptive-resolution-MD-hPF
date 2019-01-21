import numpy as np
from numpy import floor
from scipy.interpolate import RegularGridInterpolator

def test_example_field_periodic() :
    # Generate random positions.
    N = 20
    box = np.array([5.0, 5.0, 5.0])
    p = np.random.uniform(low=0.0, high=box[0], size=(N,3))
    
    # Field parameters
    m           = 3
    l           = box[0] / m
    cell_volume = l * l * l
    x           = np.linspace(0,box[0]-l,m)
    density     = np.zeros(shape=(m,m,m))

    corners = [[0,0,0], 
               [1,0,0],
               [0,1,0],
               [0,0,1],
               [1,1,0],
               [1,0,1],
               [0,1,1],
               [1,1,1]]
    corners  = np.array(corners)
    opposite = np.arange(7,-1,-1)
    
    for i in range(N) :
        x_cell = int(floor(p[i,0] / l))
        y_cell = int(floor(p[i,1] / l))
        z_cell = int(floor(p[i,2] / l))

        x = p[i,0] - x_cell*l
        y = p[i,1] - y_cell*l
        z = p[i,2] - z_cell*l

        partial_volume = np.zeros(8)
        density_ind      = [x_cell,   y_cell,   z_cell]
        density_ind_next = [x_cell+1, y_cell+1, z_cell+1]

        for j in range(8) :
            c_x = corners[j,0]
            c_y = corners[j,1]
            c_z = corners[j,2]

            partial_volume[j] = (  abs((c_x * l) - x) 
                                 * abs((c_y * l) - y)
                                 * abs((c_z * l) - z) ) / cell_volume

        if x_cell == m-1 :
            density_ind_next[0] = 0
        if y_cell == m-1 :
            density_ind_next[1] = 0
        if z_cell == m-1 :
            density_ind_next[2] = 0

        d   = density_ind
        d_n = density_ind_next

        density[d_n[0], d_n[1], d_n[2]] += partial_volume[0] # 0,0,0
        density[d  [0], d_n[1], d_n[2]] += partial_volume[1] # 1,0,0
        density[d_n[0], d  [1], d_n[2]] += partial_volume[2] # 0,1,0
        density[d_n[0], d_n[1], d  [2]] += partial_volume[3] # 0,0,1
        density[d  [0], d  [1], d_n[2]] += partial_volume[4] # 1,1,0
        density[d  [0], d_n[1], d  [2]] += partial_volume[5] # 1,0,1
        density[d_n[0], d  [1], d  [2]] += partial_volume[6] # 0,1,1
        density[d  [0], d  [1], d  [2]] += partial_volume[7] # 1,1,1

    print_fort5(box, p)
    print_adap_pos(p)

    for i in range(m) :
        for j in range(m) :
            for k in range(m) :
                print(i+1, j+1, k+1, end="")
                print("%8.2f %8.2f %8.2f %20.15f" % (i*l, j*l, k*l, density[i,j,k]))


def print_fort5(box_size, coordinates) :
    # Parameters
    coordinates          = coordinates.T
    n_molecules          = coordinates.shape[1]
    n_atoms_per_molecule = 1
    atom_label           = 'A'
    atom_id              = 1
    n_bonds              = 0
    max_n_bonds          = 6

    print("# Simulation box: x, y, z, binc")
    for i in range(3) :
        print(box_size[i],end=" ")
    print("0.0")

    print("# Total number of molecules: N")
    print(n_molecules)

    #coordinates = np.random.uniform(low=0.0, high=1.0, size=(3, n_molecules))

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


def print_adap_pos(positions) :
    N = positions.shape[0]
    for i in range(N) :
        print("positions(:,",i+1,") = [", end="")
        print(positions[i,0], end="_real64,")
        print(positions[i,1], end="_real64,")
        print(positions[i,2], end="_real64]\n")





if __name__ == '__main__':
    np.random.seed(9185)
    test_example_field_periodic()