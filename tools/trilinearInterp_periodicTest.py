import numpy as np
from numpy import floor
from scipy.interpolate import RegularGridInterpolator

def test_example_field_periodic() :
    # Generate random positions.
    N = 10
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
        density[d  [0], d_n[1], d_n[2]] += partial_volume[0] # 1,0,0
        density[d_n[0], d  [1], d_n[2]] += partial_volume[0] # 0,1,0
        density[d_n[0], d_n[1], d  [2]] += partial_volume[0] # 0,0,1
        density[d  [0], d  [1], d_n[2]] += partial_volume[0] # 1,1,0
        density[d  [0], d_n[1], d  [2]] += partial_volume[0] # 1,0,1
        density[d_n[0], d  [1], d  [2]] += partial_volume[0] # 0,1,1
        density[d  [0], d  [1], d  [2]] += partial_volume[0] # 1,1,1

    





if __name__ == '__main__':
    test_example_field_periodic()