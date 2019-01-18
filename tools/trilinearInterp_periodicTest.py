import numpy as np
from scipy.interpolate import RegularGridInterpolator

def test_example_field_periodic() :
    # Generate random positions.
    N = 10
    p = np.random.uniform(low=0.8, high=1.0, size=(N,3))
    x = np.array([0.8, 1.0])
    corners = [[0,0,0], 
               [1,0,0],
               [0,1,0],
               [0,0,1],
               [1,1,0],
               [1,0,1],
               [0,1,1],
               [1,1,1]]
    corners  = np.array(corners)
    vertices = 0.8 + 0.2*corners
    opposite = np.arange(7,-1,-1)
    density  = np.zeros(shape=(2,2,2))
    l = 0.2
    
    for i in range(N) :
        x = p[i,0]
        y = p[i,1]
        z = p[i,2]

        partial_volume    = np.zeros(8)
        for j in range(8) :
            partial_volume[j] = ((corners[j,0] - x) * 
                                 (corners[j,1] - y) * 
                                 (corners[j,2] - z))
            density[corners[opposite[j]]] += partial_volume[j]

    print(density)
        




if __name__ == '__main__':
    test_example_field_periodic()