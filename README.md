# iricore

A Python interface to IRI-2016 and IRI-2020 using `ctypes` communication.

**The extensive documentation is available at the [RTD wbsite](https://iricore.readthedocs.io/en/latest).**

## Installation

This package proved to work under Linux only (due to compilation difficulties in Windows).
If you are using Windows - consider installing [WSL](https://docs.microsoft.com/en-us/windows/wsl/install).

### Prerequisites

- CMAKE

```
sudo apt instal cmake
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