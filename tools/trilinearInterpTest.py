import numpy as np
from scipy.interpolate import RegularGridInterpolator

# Generate random field values.
field  = np.random.uniform(size=(2,2,2)) 
print("\nField values:")
for i in range(2) :
    for j in range(2) :
        for k in range(2) :
            print("density_field(", i+1, ",", j+1, ",", k+1, ") = ", end="")
            print(field[i,j,k], end="")
            print("_real64")

# Scipy interpolator.
x = np.array([0, 1])
interpolator = RegularGridInterpolator((x,x,x), field)

# Random points to evaluate at.
N = 10
points = np.zeros(shape=(N, 3))
for i in range(N) :
    points[i] = np.random.uniform(size=3)

print("\nEvaluation points:")
for i in range(N) :
    print("points(:,", i+1, ") = [", end="")
    for j in range(3) :
        print(points[i,j], end="")
        if j != 2 :
            print("_real64, ", end="")
        else :
            print("_real64]")


# Interpolated values.
v = interpolator(points, method="linear")

print("\nInterpolated values:")
for i in range(N) :
    print("values(", i+1, ") = ", end="")
    print(v[i], end="")
    print("_real64")

