# mfpatbc
### Version 2.0.0

## Introduction
mfpatbc is a Python module to create input data for a phase-averaged tidal boundary (PA-TBC) condition for MODFLOW models created with [flopy](https://github.com/modflowpy/flopy). The PA-TBC is explained in detail in **Haehnel et al. (202X)**.

Supported versions of MODFLOW are MODFLOW-2005, MODFLOW-NWT, SEAWAT, and MODFLOW 6 (constant density). The module defines stress period based input data for the General Head Boundary (GHB) and the Drainage Boundary (Package) as required for the PA-TBC. For MODFLOW-2005, defining freshwater heads for the GHB package as required by the Seawater Intrusion Package (SWI2) is supported.

For use with SEAWAT, the GHB input for the Source and Sink Mixing (SSM) package of MT3DMS can be defined.

Should you find any error or have ideas for improvements, I would be very grateful if you let me know in 'Issues' or 'Discussions'.

## Installation
### Python module 
Please download the source code.
Place it: 
* in a directory where Python looks for modules (sys.path)
* local in the directory of the script that loads the module
* in any directory which is added to the sys.path in your script 
```python
>>> sys.path.append('path/to/your/directory')
```
Or any other, more elegant method you may know. 

### Compiling Fortran code
The module has a functionality to write GHB and DRN input files independent of flopy (mfpatbc.write_ghb() and mfpatbc.write_drn()) using some Fortran code to write the stress period data. The looping through the stress periods while writing input files is still in Python but the writing of the cell information per stress period is outsourced to Fortran. This is especially time-saving for models with many PA-TBC cells and/or many stress periods. Depending on the model the runtime savings can be around 50 % compared to the flopy writing utilities for the respective input files.

The connection between Python and Fortran is obtained by numpy.f2py ([documentation](https://numpy.org/doc/stable/f2py/index.html)). The Fortran script "mfpatbc/write_spd_fort.f90" needs to be compiled on the users operating system in order for the file writing capabilities of mfpatbc to work. This needs some work on the console or terminal (Considering you probably want this, because you run a quite large model, you're likely in the game, so I guess this is fine ;) ).

Steps:
1. Open terminal or console of your operating system
2. Navigate to the folder 'mfpatbc/' that contains the scripts of this module
3. Compile the Fortran script by entering (if there are any problems on Windows there is some information on numpy.f2py with windows in the [documentation](https://numpy.org/doc/stable/f2py/windows/index.html))
```bas
$ python -m numpy.f2py -c write_spd_fort.f90 -m write_spd_fort
```
The resulting file provides an interface between Python and Fortan so that the subroutine in the Fortran script, that does the heavy lifting, can be accessed just as a normal Python module by mfpatbc. All this is only tested on an Ubuntu system, so errors may (or should?) occure elsewhere.


## Documentation
Is provided within the methods, functions, and classes as docstring. Access e.g. via
```python
help(mfpatbc.PATBC)
?mfpatbc.PATBC
??mfpatbc.PATBC
```

## Getting started
The folder 'examples' provides Jupyter Notebooks showcasing the use of the module for MODFLOW-2005, SEAWAT, and MODFLOW 6 with flopy.

## How to cite
Haehnel, Patrick (2023). mfpatbc v2.0.0 [Software] URL: https://github.com/phaehnel/mfpatbc

## References
**Haehnel et al. (202X)**
