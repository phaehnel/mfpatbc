#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 13:10:21 2023

@author: patrick
"""

import os
import numpy as np
import mfpatbc
from .write_spd_fort import *

def define_filename(m, filename, extension):
    """
    Set filename for boundary condition input files. Filename input not provided
    model name is used. If filename is provided it is checked for presence of
    file extension. Is added if not present.

    Parameters
    ----------
    m : flopy model object
        Either MODFLOW-2005, SEAWAT or MODFLOW6 model class object.
    filename : str
        Intended name of the file.
    extension : str
        Intended file extension.

    Returns
    -------
    filename : str
        Filename for boundary condition input file.

    """
    
    # If filename is not defined use name of model object and append file 
    # extension
    if filename is None:
        filename = f'{m.name}.{extension}'
    # Otherwise set user defined filename and check if file extension is
    # present
    else:
        has_extension = extension in filename
        filename = filename if has_extension else f'{filename}.{extension}'
    
    return filename


def get_maximum_number_bc_cells(stress_period_data):
    """
    Find maximum number of boundary condition cells for all stress periods.

    Parameters
    ----------
    stress_period_data : dictionary
        Dictionary with stress period data of any kind.

    Returns
    -------
    max_ncells : int
        Maximum number of boundary condition cells.

    """
    
    max_ncells = 0
    # Check number of entries in each stress period and update maximum number
    # found so far if necessary.
    for values in stress_period_data.values():
        ncells = len(values)
        max_ncells = max_ncells if ncells <= max_ncells else ncells
        
    return max_ncells


def get_options_string(options, precision = None):
    """
    Create the options string for auxiliary variables for MODFLOW boundary
    condition input files.

    Parameters
    ----------
    options : str or list of strs
        Names of auxiliar variables to write to MODFLOW boundary condition
        input file.
    precision : int, optional
        length of each string in the return string (default is None). Only
        required if options is list of strings with len(options) > 1.

    Returns
    -------
    options_str : str
        String of auxiliar variables to write to MODFLOW boundary condition
        input file.

    """
    
    is_string = isinstance(options, str)
    
    # If already string do nothing
    if is_string:
        options_str = f' {options:{precision}s}'
        
    # If one element list return entry
    elif not is_string and (len(options) == 1):
        options_str = f' {options[0]:{precision}s}'
    
    # If multi element list create string with elements separated by delimiter.
    else:
        options_str = ''
        for value in options:
            options_str += f' {value:{precision}s}'
    
    return options_str
    

def write_drn(m, stress_period_data, mxactb = None, idrncb = 0, options = None,
              filename = None, filepath = None, verbose = False, 
              precision = {'cell': 9, 'data': 15}, newline = '\n'):
    """
    Write Drainage Boundary Package (DRN) input file based on stress period data 
    defined e.g. with method PATBC.get_drn_stress_period_data(). 
    
    This function allows to write auxiliar variables to the stress period data
    input of the DRN input file.
    
    See section 'Examples' in the documentation of mfpatbc.PATBC() for the 
    concept on how to use this function.
    
    So far, only supports writing of stress_period_data dictionaries in 
    format required for flopy.modflow.ModflowDrn (MODFLOW-2005 or SEAWAT).
    
    This function relies on a Fortran script to write the stress period data. 
    Cf. to the documentation of the mfpatbc module for 
    information how to create the Fortran executable required. 
    Check function write_stress_period_data for details how stress period data
    is written to the input file.
    
    Input data must meet requirements defined for DRN input in Harbaugh (2005).

    Parameters
    ----------
    m : flopy model object
        Either MODFLOW-2005 or SEAWAT model class object.
    stress_period_data : dict
        Stress period data for DRN. Format of dictionary entries is 
        (layer, row, column, stage, conductance, AUX), where AUX is any
        number of auxiliary variables defined in attribute options.
    mxactb : int, optional
        Maximum number of DRN cells required in simulation. 
        If not provided number is determined from stress period data by function
        get_maxmimum_number_bc_cells. The default is None.
    idrncb : int, optional
        Flag determining if cell-by-cell flow term will be written to output. 
        The default is 0.
    options : str or list of strs, optional
        Names of auxiliar variables to write to MODFLOW boundary condition
        input file. Either provide list of strings, which are then formated
        accoring to precision['cell'], or provide string which can be written as
        is to DRN input file. The default is None.
    filename : str, optional
        Name of the DRN input file. If not provided, modelname of the flopy
        model object is used. The default is None.
    filepath : TYPE, optional
        Directory to write DRN input file. If not provided, model workspace of
        the flopy model object is used. If defined make sure that the directory 
        is defined accordingly in the name file. The default is None.
    verbose : bool, optional
        If True, writes information on progress to console. The default is False.
    precision : dict, optional
        precision for cell and data information, when writing DRN input file. 
        'cell' refers to the layer, row, and column information, 
        'data' refers to stage, conductivity and auxiliary variables.
        So far, this only effects lines 0 and 2 of the input file.
        The stress period data precision is defined in the fortran script itself.
        The default is {'cell': 9, 'data': 15}.
    newline : str, optional
        Newline command. Should typically not be changed.
        
    Returns
    -------
    None.
    
    References
    ----------
    Harbaugh, A. W. (2005). MODFLOW-2005, The U.S. Geological Survey Modular 
    Ground-Water Model—the Ground-Water Flow Process (Techniques and Methods 
    No. 6-A16). Reston: U.S. Geological Survey.
    
    """
    
    # Define file name
    filename = define_filename(m, filename, 'drn')
    
    # Set file location
    if filepath is None:
        filepath = m.model_ws
    
    # Set maximum number of GHB boundary cells
    if mxactb is None:
        mxactb = get_maximum_number_bc_cells(stress_period_data)
        
    # Define options string
    options_str = '' if options is None else\
        get_options_string(options, precision['cell'])
        
    if verbose:
        print('Start writing DRN input file')
        
    # Open file
    fname = os.path.join(filepath, filename)
    f = open(fname, 'w')
    
    # Line numbering according to Harbaugh (2005)
    # Line 0
    f.write(f'# DRN package, generated by mfpatbc version {mfpatbc.__version__}' +
            f'{newline}')
    
    # Line 2
    f.write(
        f' {mxactb:{precision["cell"]}d} {idrncb:{precision["cell"]}d}' +\
        f'{options_str}{newline}'
    )
    
    f.close()
        
    # Lines 5 and 6
    write_stress_period_data(fname, m.nper, stress_period_data, verbose)
    
    if verbose:
        print('Finished writing DRN input file\n')
 
      
def write_ghb(m, stress_period_data, mxactb = None, ighbcb = 0, options = None,
              filename = None, filepath = None, verbose = False, 
              precision = {'cell': 9, 'data': 15}, newline = '\n'):
    """
    Write General Head Boundary Package (GHB) input file based on stress period 
    data defined e.g. with method PATBC.get_ghb_stress_period_data(). 
    
    This function allows to write auxiliar variables to the stress period data
    input of the GHB input file.
    
    See section 'Examples' in the documentation of mfpatbc.PATBC() for the 
    concept on how to use this function.
    
    So far, only supports writing of stress_period_data dictionaries in 
    format required for flopy.modflow.ModflowGhb (MODFLOW-2005 or SEAWAT).
    
    This function relies on a Fortran script to write the stress period data. 
    Cf. to the documentation of the mfpatbc module for 
    information how to create the Fortran executable required. 
    Check function write_stress_period_data for details how stress period data
    is written to the input file.
    
    Input data must meet requirements defined for GHB input in Harbaugh (2005).

    Parameters
    ----------
    m : flopy model object
        Either MODFLOW-2005 or SEAWAT model class object.
    stress_period_data : dict
        Stress period data for GHB. Format of dictionary entries is 
        (layer, row, column, stage, conductance, AUX), where AUX is any
        number of auxiliary variables defined in attribute options.
    mxactb : int, optional
        Maximum number of GHB cells required in simulation. 
        If not provided number is determined from stress period data by function
        get_maxmimum_number_bc_cells. The default is None.
    ighbcb : int, optional
        Flag determining if cell-by-cell flow term will be written to output. 
        The default is 0.
    options : str or list of strs, optional
        Names of auxiliar variables to write to MODFLOW boundary condition
        input file. Either provide list of strings, which are then formated
        accoring to precision['cell'], or provide string which can be written as
        is to GHB input file. The default is None.
    filename : str, optional
        Name of the GHB input file. If not provided, modelname of the flopy
        model object is used. The default is None.
    filepath : TYPE, optional
        Directory to write GHB input file. If not provided, model workspace of
        the flopy model object is used. If defined make sure that the directory 
        is defined accordingly in the name file. The default is None.
    verbose : bool, optional
        If True, writes information on progress to console. The default is False.
    precision : dict, optional
        precision for cell and data information, when writing GHB input file. 
        'cell' refers to the layer, row, and column information, 
        'data' refers to stage, conductivity and auxiliary variables. 
        So far, this only effects lines 0 and 2 of the input file.
        The stress period data precision is defined in the fortran script itself.
        The default is {'cell': 9, 'data': 15}.
    newline : str, optional
        Newline command. Should typically not be changed.
        
    Returns
    -------
    None.
    
    References
    ----------
    Harbaugh, A. W. (2005). MODFLOW-2005, The U.S. Geological Survey Modular 
    Ground-Water Model—the Ground-Water Flow Process (Techniques and Methods 
    No. 6-A16). Reston: U.S. Geological Survey.
    
    """
    
    # Define file name
    filename = define_filename(m, filename, 'ghb')
    
    # Set file location
    if filepath is None:
        filepath = m.model_ws
    
    # Set maximum number of GHB boundary cells
    if mxactb is None:
        mxactb = get_maximum_number_bc_cells(stress_period_data)
        
    # Define options string
    options_str = '' if options is None else\
        get_options_string(options, precision['cell'])
        
    if verbose:
        print('Start writing GHB input file')
        
    # Open file
    fname = os.path.join(filepath, filename)
    f = open(fname, 'w')
    
    # Line numbering according to Harbaugh (2005)
    # Line 0
    f.write(f'# GHB package, generated by mfpatbc version {mfpatbc.__version__}' +
            f'{newline}')
    
    # Line 2
    f.write(
        f' {mxactb:{precision["cell"]}d} {ighbcb:{precision["cell"]}d}' +\
        f'{options_str}{newline}'
    )

    f.close()
        
    # Lines 5 and 6
    write_stress_period_data(fname, m.nper, stress_period_data, verbose)
        
    if verbose:
        print('Finished writing GHB input file\n')
    
    
def write_stress_period_data(fname, nper, stress_period_data, verbose):
    """
    Write stress period data in format of stress_period_data dictonaries of
    flopy.modflow stress period based boundary condition packages 
    (MODFLOW-2005 or SEAWAT).
    
    The writing of the stress period data (Lines 5 and 6 of input file) relies
    on a fortran script which is faster than the flopy implementation of the
    file writing operation for GHB and DRN. The Fortran script is wrapped for
    use in Python using numpy.f2py, which allows to access the Fortran code
    like a Python module. Cf. to the documentation of the mfpatbc module for 
    information how to create the Fortran executable required.
    
    No data is written for stress periods with smaller numbers than first key 
    in stress_period_data dictionary. When stress period number is larger than
    first key but no entry is in the stress_period_data dictionary, data from
    the last stress period defined in the stress_period_data dictionary
    is reused.

    Parameters
    ----------
    f : file object
        File object of open() to write data into.
    nper : int
        number of stress periods.
    stress_period_data : dict
        Stress period data for boundary condition package. Format of dictionary 
        entries is (layer, row, column, stage, conductance, AUX), where AUX is any
        number of auxiliary variables defined by the calling function.
    verbose : bool
        If True, writes information on progress to console.

    Returns
    -------
    None.

    """
    
    
    keys = list(stress_period_data.keys())
    first_per = keys[0]
    for kper in range(nper):
        
        if verbose:
            print(f'Write stress period {kper}')
        
        # Stress periods before first dictionary entries or empty dict entries
        if (kper < first_per):
            itmp = 0
            do_write = False
            spd = np.asfortranarray([[-9999, -9999, -9999, -9999, -9999]])
            
        # If kper present in keys write data to file
        elif kper in keys:
            spd_kper = stress_period_data[kper]
            itmp = len(spd_kper)
            
            spd = np.asfortranarray(spd_kper.tolist())#.astype(np.float32)
            # Add 1 to indices of lay, row, col
            if itmp > 0:
                spd[::, :3] += 1
                do_write = True
            else:
                do_write = False
                spd = np.asfortranarray([[-9999, -9999, -9999, -9999, -9999]])
            
        # If kper not in dictionary keys, reuse previous stress period data
        else:
            itmp = -1
            do_write = False
            spd = np.asfortranarray([[-9999, -9999, -9999, -9999, -9999]])
            
        nrow, ncol = spd.shape

        write_spd_fort(
            fname, spd, itmp, kper + 1, do_write, nrow, ncol
        )
