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

