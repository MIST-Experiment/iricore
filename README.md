# iricore
A fork of [iri2016](https://github.com/space-physics/iri2016). `iricore` implements a couple of optimizations to the `iri2016`
core Fortran code and Python interface to make it faster:
1. Optimization of data files reading gives boost in case of iteration over large list of coordinates;
2. `f2py` interface provides faster communication between Python and Fortran.

Overall, this gives performance boost up to ~100x (see `examples/comparison.py`).

**Important!** Because this package is mainly used for the [MIST experiment](http://www.physics.mcgill.ca/mist/), 
the `iricore` cuts off calculation of unnecessary atmospheric parameters available in `iri2016`, leaving only electron density
and electron temperature. All other parameters can be restored on demand (please contact me).

## Installation

This package proved to work under Linux only (due to compilation difficulties in Windows). 
If you are using Windows - consider installing [WSL](https://docs.microsoft.com/en-us/windows/wsl/install).

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
Now you can simply install it via `pip`:

```
python3 -m pip install iricore
```

## Data files
`IRI2016` model depends on [data files](http://irimodel.org/indices/) which are regularly updated.
`iricore` does not autoupdate those, but provides tool for quick update. You can run from terminal
```
python3 -c "import iricore; iricore.update()"
```

or add

```
import iricore
iricore.update()
```
to any Python script.

## Usage
For usage examples see `examples/`.