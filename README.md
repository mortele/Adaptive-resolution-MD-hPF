[![Build Status](https://travis-ci.com/mortele/Adaptive-resolution-MD-hPF.svg?token=81VUNKkUYjZSicZzs1NR&branch=master)](https://travis-ci.com/mortele/Adaptive-resolution-MD-hPF) 
[![codecov](https://codecov.io/gh/mortele/Adaptive-resolution-MD-hPF/branch/master/graph/badge.svg?token=ayq0rwnrot)](https://codecov.io/gh/mortele/Adaptive-resolution-MD-hPF) 
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

## Adaptive resolution MD-hPF
Hybrid particle-field method based on [the 2009 paper by Milano & Kawakatsu](https://aip.scitation.org/doi/abs/10.1063/1.3142103). The adaptive resolution interface is inspired by [the 2005 paper by Praprotnik, Delle Site, & Kremer](https://aip.scitation.org/doi/pdf/10.1063/1.2132286). The molecular dynamics is written loosely following [the book of Frenkel & Smith](https://www.elsevier.com/books/understanding-molecular-simulation/frenkel/978-0-12-267351-1).

### Dependencies
Building the binaries from source requires 
* the GNU Fortran compiler (GFortran) *v ≥ 4.8.4*
* GNU make *version ≥ 3.81*

### How to build
Clone the repository,
```
> git clone git@github.com:mortele/Adaptive-resolution-MD-hPF.git MDhPF
```
Navigate into the directory and run `make` to compile, and/or `make test` to compile the test suite,
```
> cd MDhPF
> make
> ./bin/AdapResoMD-hPF.app
```
or
```
> cd MDhPF
> make test
> ./bin/UnitTests.app
```

### TODO
* Finish testing `file_writer.f90`
    - Semi-automate visualization with VMD(?)
* Test `compute_density_field()` to make sure a variable number of field nodes in each direction works properly
* Finish testing `compute_density_gradient()`
* Perform **larger than unit-test scale** testing of MD code and FCC lattice
* Implement cell lists and neighbor lists
* Implement *hybrid particle-field* force calculation and test energies
* Implement higher order central difference schemes for calculating the gradient density
* Test with very high density field resolution that we recreate all-atom MD dynamics (is higher order derivative schemes neccessary for this?)
* Experiment with how high field resolution is needed before FCC lattice is energetically stable in the hPF case
* Implement adaptive resolution hPF-MD interface 
    - Test energy conservation as particles diffuse across the boundary
    - Test that net drift across the boundary in both directions is zero

