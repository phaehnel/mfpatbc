# mfpatbc
### Version 1.0.0

## Introduction
mfpatbc is a Python module to create input data for a phase-averaged tidal boundary (PA-TBC) condition for MODFLOW models created with [flopy](https://github.com/modflowpy/flopy). The PA-TBC is explained in detail in **Haehnel et al. (202X)**.

Supported versions of MODFLOW are MODFLOW-2005, MODFLOW-NWT, SEAWAT, and MODFLOW 6 (constant density). The module defines stress period based input data for the General Head Boundary (GHB) and the Drainage Boundary (Package) as required for the PA-TBC.

For use with SEAWAT, the GHB input for the Source and Sink Mixing (SSM) package of MT3DMS can be defined. The modul also provides functions to write GHB and DRN input files with SEAWAT auxiliary variables IGHBELEV and IDRNELEV, which is not supported by flopy yet.

Should you find any error or have ideas for improvements, I would be very grateful if you let me know in 'Issues' or 'Discussions'.

## Installation

Please download the source code.
Place it: 
* in a directory where Python looks for modules (sys.path)
* local in the directory of the script that loads the module
* in any directory which is added to the sys.path in your script 
```python
>>> sys.path.append('path/to/your/directory')
```
Or any other, more elegant method you may know. 

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
Haehnel, Patrick (2020). mfpatbc v1.0.0 [Software] URL: https://github.com/phaehnel/mfpatbc

## References
**Haehnel et al. (202X)**