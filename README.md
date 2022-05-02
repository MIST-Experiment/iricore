# iricore
A fork of [iri2016](https://github.com/space-physics/iri2016). `iricore` implements a couple of optimizations to the `iri2016`
core Fortran code and Python interface to make it faster:
1. Optimization of data files reading gives boost in case of iteration over large list of coordinates;
2. `f2py` interface provides faster communication between Python and Fortran.

Overall, this gives up to ~100x performance boost (see `examples/comparison.py`).

**Important!** Because this package is mainly used for the [MIST experiment](http://www.physics.mcgill.ca/mist/), 
the `iricore` cuts off calculation of unnecessary atmospheric parameters available in `iri2016`, leaving only electron density
and electron temperature. All other parameters can be restored on demand (please contact me).

## Installation

This package proved to work under Linux only (due to compilation difficulties in Windows). 
If you are using Windows - consider isntalling [WSL](https://docs.microsoft.com/en-us/windows/wsl/install).

### Prerequisites
- Git
```
sudo apt instal git
```

- Fortran compiler, e.g. `gfortran`
```
sudo apt isntall gfortran
```

### Installing package
[//]: # ()
[//]: # (Then you can install the package via pip)

```
python3 -m pip install git+https://github.com/lap1dem/iricore
```

## Usage
For usage examples see `examples/`.