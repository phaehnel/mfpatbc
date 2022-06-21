#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Patrick Haehnel
2022-05-18
"""

import os
import numpy as np
import pandas as pd
import flopy
import mfpatbc


class PATBC():
    """
    Object to create stress period data input for MODFLOW General Head (GHB) and
    Drainage (DRN) boundary package in flopy. Object holds and creates input 
    data for one single stress period. See 'Notes' for background information 
    and details. See 'Examples' for concept of implementation with flopy.
    
    Methods of this class can be used to retrieve input for dictonary 
    of attribute 'stress_period_data' in flopy objects for GHB and DRN 
    boundary conditions.
    
    The module allows to create respective output for MODFLOW-2005, SEAWAT and
    MODFLOW6. Implementation is only done for structured grid in MODFLOW6.
    For MODFLOW-2005, the module is also able to caluculate boundary heads as 
    required by the Seawater Intrusion Package (SWI2) (Bakker et al., 2013). See
    attribute version.
    
    Hand all parameters in unit of the MODFLOW model to ensure correct 
    calculation of vertical conductance for the GHB and DRN.
    
    Attributes (user defined)
    -----
    m : flopy model object
        Either MODFLOW-2005, SEAWAT, or MODFLOW6 model class object.
    hs : float or array of floats (nrow, ncol)
        Phase-averaged sea level.
    A : float or array of floats (nrow, ncol)
        Tidal amplitude.
    T : float or array of floats (nrow, ncol)
        Period length of tidal cycle.
    hk : float or array of floats (nrow, ncol)
        Horizontal hydraulic conductivity, dimensional [m/s].
    vka : float or array of floats (nrow, ncol)
        vertical anisotropy.
    sy : float or array of floats (nrow, ncol)
        Specific yield.
    z_D : float or array of floats (nrow, ncol)
        Elevation of the aquifer base (e.g. m asl). Used to calculate
        saturated aquifer thickness D = hs - z_D
    slope : float or array of floats (nrow, ncol)
        beach slope as gradient, tan(beta)
    idx_top_layer : array of ints (nrow, ncol), optional
        Indices of layer whose cell top represent surface elevation 
        (default is None). If not provided data is retrieved from IBOUND 
        (MODFLOW-2005, SEAWAT) or IDOMAIN (MODFLOW6) by searching for uppermost
        active model cells (method get_top_layer). Values are used in methods
        get_drn_stress_period_data and get_ghb_stress_period_data to define
        boundary condition layer indices. Row and column combinations where
        all cells are inactive should hold a index < 0. This ensures correct
        set up of the attribute active.
    allow : array of bools (nrow, ncol), optional
        True if cell is allowed to hold a phase-averaged tidal boundary condition
        (default is None). If not provided all cells will be allowed to hold
        a phase-averaged tidal boundary condition and each cell is assigned either
        a GHB or a DRN.
    version : str, optional
        Define version of of MODFLOW (default is None). Should only be defined
        by user if MODFLOW-2005 is used with SWI2 package (seawater intrusion)
        package. Then set to 'swi2'. This allows calculation of freshwater heads
        as required by SWI2.
    rho_fresh : float, optional
        Density of fresh water (default is 1000 kg/m³). Required for SEAWAT
    rho_salt : float, optional
        Density of salt water (default is 1025 kg/m³). Required for SEAWAT
    
    
    Attributes (derived from flopy model object)
    -----
    version : str
        MODFLOW version. Required for correct access of certain attributes of
        flopy model object
    nlay : int
        number of model layers, derived from m.dis.botm
    nrow : int
        number of model rows, m.dis.nrow
    ncol : int
        number of model columns, m.dis.ncol
    delr : array of floats (ncol)
        model cell spacings along rows, m.dis.delr.array 
    delc : array of floats (nrow)
        model cell spacings along columns, m.dis.delc.array
    
    
    Attributes (derived from methods)
    -----
    active : array of bools (nrow, ncol)
        True if cell location holds a phase-averaged boundary condition in any
        model layer. Determined based on setting of idx_top_layer
    z : array of floats (nrow, ncol)
        Surface elevations, derived from m.dis.top based on method get_top_layer.
        Defined as upper cell boundary of first active layer in (nrow, ncol)
        if idx_top_layer = None. Ideally, this value should coincide with surface
        elevations e.g. from a digital elevation model.
    bottom : array of floats (nrow, ncol)
        Bottom elevation of layer where cell with surface elevation at cell top
        is found. Determined based on setting of idx_top_layer
    vertical_conductance : array of floats (nrow, ncol)
        Vertical conductance (C_v = vk*A/l) of the phase-averaged boundary 
        condition. Derived by method get_vertical_conductance
    lay_top : array of ints
        layer of cells holding PA-TBC
    row_top : array of ints
        row of cells holding PA-TBC
    col_top : array of ints
        column of cells holding PA-TBC
    active : array of bools (nrow, ncol)
        True, if any layer holds PA-TBC
    
    
    Methods
    -----
    get_drn_stress_period_data :
        Hands back numpy.array required as input for the stress_period_data
        dictionary of the DRN package for the given stress period.
    get_drn_stress_period_data :
        Hands back numpy.array required as input for the stress_period_data
        dictionary of the GHB package for the given stress period.
    get_freshwater_heads :
        Hands back freshwater heads calculated according to equation 6 in 
        Post et al. (2007). Assumes GHB reservoir with density rho_salt and
        bottom of reservior at surface elevation (aquifer-ocean interface).
    get_ssm_stress_period_data :
        Hands back numpy.array required as input for the stress_period_data
        dictionary of the SSM package to define concentrations within the GHB
        reservoirs of the PA-TBC.
    get_surface_elevation :
        Hands back surface elevation of the aquifer ocean interface based on
        definition of attribute idx_top_layer.
    get_top_layer :
        Hands back numpy.array holding information which layer holds the cell
        whose top boundary represents surface elevation. Is determined based on
        IBOUND (MODFLOW-2005, SEAWAT) or IDOMAIN (MODFLOW6).
    get_vertical_conductance :
        Calculates vertical conductance for the phase-averaged boundary condition.
    make_array :
        Checks if attribute is scalar. If so, creates array of size (nrow, ncol)
        for this attribute.
    
    
    Returns
    -----
    Object of class PATBC
    
    Examples
    -----
    **Constant density simulations**
    
    Setup flopy model with all required packages.
    
    >>> m = flopy.modflow.Modflow(...)
    
    Loop through stress periods. Create PATBC object for every stress period
    and write stress period data for GHB and DRN to dictionary. The NOPRINT
    option is very convinient for large models to prevent writing all
    GHB and DRN info into the list file.
    
    >>> ghbspd = {}
    >>> drnspd = {}
    >>> for k in range(m.nper):
    >>>     patbc = mfpatbc.PATBC(...)
    >>>     ghbspd[k] = patbc.get_ghb_stress_period_data()
    >>>     drnspd[k] = patbc.get_drn_stress_period_data()
    >>> ghb = flopy.modflow.ModflowGhb(m, stress_period_data = ghbspd, 
                                       option = ['NOPRINT'])
    >>> drn = flopy.modflow.ModflowDrn(m, stress_period_data = drnspd, 
                                       option = ['NOPRINT'])


    **Sharp-interface simulations with SWI2 package**
    
    Same as for constant density simulations. Only requires definition of 
    attribute version = 'swi2' in mfpatbc.PATBC() object
    which enables calculation of freshwater heads for the GHB.
    
    
    **Variable density simulations**
    
    Setup flopy model with all required packages.
    
    >>> swt = flopy.seawat.Seawat(...)
    
    Loop through stress periods. Create PATBC object for every stress period
    and write stress period data for GHB and DRN to dictionary. Create dummy
    objects for GHB and DRN in flopy, as auxiliary variables IGHBELEV and
    IDRNELEV for SEAWAT are not supported by flopy. Dummy objects will make
    flopy write GHB and DRN to name file. Boundary heads of GHB are calculated
    as freshwater heads, assuming the density set for the GHB reservior.
    
    >>> ghbspd = {}
    >>> drnspd = {}
    >>> for k in range(swt.nper):
    >>>     patbc = mfpatbc.PATBC(...)
    >>>     ghbspd[k] = patbc.get_ghb_stress_period_data()
    >>>     drnspd[k] = patbc.get_drn_stress_period_data()
    >>> ghb = flopy.modflow.ModflowGhb(swt, stress_period_data = {0: [0, 0, 0, 0, 0]})
    >>> drn = flopy.modflow.ModflowDrn(swt, stress_period_data = {0: [0, 0, 0, 0, 0]})
    
    Create all MT3DMS packages and SEAWAT VDF package.
    For Source and Sink Mixing package (SSM) define GHB concentrations.
    Can use the last patbc object created as SSM dictionary for GHB is created
    independent on stress period data. Sets SSM GHB to all cells at aquifer-ocean
    interface (surface elevation) that are active and allowed to hold a PA-TBC.
    Define SSM stress period dependent if attribute allow changes during simulation.
    Make sure that maximum number of SSM cells (mxss) is set large enough in 
    flopy SSM object.
    
    >>> ssmghbspd = patbc.get_ssm_stress_period_data(conc_ghb = 35)
    >>> ssm = flopy.mt3d.Mt3dSsm(swt, stress_period_data = ssmghbspd, mxss = ...)
    
    Write input files. First remove all files in model workspace, otherwise some
    strange behavior with MT3DMS output binary files may occure. Write all 
    packages but GHB and DRN. This saves time for simulations with many stress
    periods. Write GHB and DRN input files separately to allow inclusion of 
    auxiliary variables of SEAWAT.
    
    >>> filenames = glob.glob(swt.model_ws + '/*')
    >>> for f in filenames:
    >>>     try:
    >>>         os.remove(f)
    >>>     except:
    >>>         print(f + ' could not be removed.')
    >>> packages = swt.get_package_list()
    >>> packages = [name for name in packages if (name != 'GHB') & (name != 'DRN')]
    >>> swt.write_input(SelPackList = packages)
    >>> mfpatbc.write_ghb(swt, spd_ghb, options = ['NOPRINT', 'IGHBELEV'])
    >>> mfpatbc.write_drn(swt, spd_drn, options = ['NOPRINT', 'IDRNELEV'])
    
    
    Notes
    -----
    **Phase-averaged tidal boundary condition**
    
    The phase-averaged tidal boundary condition (PA-TBC) is 
    based on an analytical solution presented by Nuttle (1991), hlim, which describes
    phase-averaged intertidal head distribution under tidal forcing for very shallow 
    sloping beaches (mfpatbc.phase_average_tbc.pavg_heads_lim).

    Empirical correction function g adapts this solutions to deviations from the
    assumptions of this solution due to parameters 
        - horizontal hydraulic conductivity (hk) 
        - specific yield (sy)
        - saturated aquifer thickness (D);
          D = hs - z_D, hs -> phase-averaged sea level, 
          z_D -> aquifer base
        - vertical anisotropy (vka = vk/hk, vk is vertical hydraulic conductivity) 
        - beach slope (beta)
    using a linear regression model describing the phase-averaged exit point 
    elevation of the groundwater table (zep). This is the elevation at which
    groundwater table is equal to surface elevation on phase-averaged basis
    (mfpatbc.phase_average_tbc.predict_zep). Note that parameters with dimensions
    (hk, D) are nondimensionalized by tidal period length T and tidal amplitude A.
    (hk*T/A and D/A). Therefore, input data does not need to be adapted to any
    specific unit as long as units are consistent between parameters.

    Adaption of hlim is a cutoff threshold where zep is used instead of hlim
    whenever hlim > zep (mfpatbc.phase_average_tbc.pavg_heads).

    All cells with surface elevations z < zep are assigned a GHB and all cells
    with z >= zep are assigned a DRN, should the respective cells be allowed
    to hold the PA-TBC. 

    Additionally, a function is included allowing a user defined, constant setting of
    the empirical correction function (mfpatbc.phase_average_tbc.pavg_heads_zep_const).

    **Variable density simulations with SEAWAT**
    
    GHB and DRN boundary conditions are defined equivalently to their definition
    in Mulligan et al. (2011) (e.g. Figure 1). Fresh water heads for boundary stage
    of GHB are retrieved by method get_freshwater_heads, but the check for cells
    being assigned a GHB is done with the point water head (which assumes 
    density rho_salt). The boundary stage of DRN is not modified compared to 
    constant density simulations.
    
    The bottom of the reservior of GHB and the bottom elevation of the drain
    (SEAWAT auxiliar variables IGHBELEV and IDRNELEV) are set to respective 
    cell values of attribute z (surface elevation, aquifer-ocean interface).
    
    Fresh water heads of the GHB are calculated with method get_freshwater_heads 
    assuming that the elevation head (z) (Post et al., 2007; eq. 6) is the
    reservoir bottom, i.e. surface elevation.
    
    **Sharp interface simulations with SWI2 package**
    
    The freshwater heads of the GHB in the PA-TBC are calculated
    at the top elevation of the uppermost active model layer for 
    each row and column combination. This definition meets the requirement of 
    the SWI2 package for head definitions in head-dependent boundary conditions
    (Bakker et al., 2013; p. 40).


    References
    -----
    Bakker, M., Schaars, F., Hughes, J. D., Langevin, C. D., 
    & Dausmann, A. M. (2013). Documentation of the Seawater Intrusion (SWI2) 
    Package for MODFLOW (Techniques and Methods No. 6-A46). 
    Reston: U.S. Geological Survey. https://doi.org/10.3133/tm6A46
    
    Mulligan, A. E., Langevin, C., & Post, V. E. A. (2011). 
    Tidal Boundary Conditions in SEAWAT. Ground Water, 49(6), 866–879. 
    https://doi.org/10.1111/j.1745-6584.2010.00788.x

    Nuttle, W. K. (1991). Comment on “Tidal dynamics of the water table in 
    beaches” by Peter Nielsen. Water Resources Research, 27(7), 1781–1782. 
    https://doi.org/10.1029/91WR00939
    
    Post, V., Kooi, H., & Simmons, C. (2007). Using Hydraulic Head Measurements 
    in Variable-Density Ground Water Flow Analyses. Ground Water, 45(6), 664–671. 
    https://doi.org/10.1111/j.1745-6584.2007.00339.x

    """
    
    def __init__(self, m, hs, A, T, hk, vka, sy, z_D, slope, idx_top_layer = None,
                 allow = None, version = None, rho_fresh = 1000, rho_salt = 1025):
        
        # Assign attributes
        self.version = m.version if version is None else version
        self.top = m.dis.top.array
        self.botm = m.dis.botm.array
        self.nlay = self.botm.shape[0]
        self.nrow = m.dis.nrow.get_data() if self.version == 'mf6' else m.dis.nrow
        self.ncol = m.dis.ncol.get_data() if self.version == 'mf6' else m.dis.ncol
        
        # If not provided all cells representing surface elevation will
        # be determined from IBOUND or IDOMAIN, respectively.
        self.idx_top_layer = idx_top_layer
        self.ibound = m.dis.idomain.array if self.version == 'mf6'\
            else m.bas6.ibound.array
        self.lay_top, self.row_top, self.col_top, self.active = self.get_top_layer()
            
        self.bottom = self.botm[self.lay_top, self.row_top, self.col_top]
        self.bottom = self.bottom.reshape((self.nrow, self.ncol))
        self.z = self.get_surface_elevation()
        self.hs = hs
        self.A = A
        self.T = T
        self.hk = hk
        self.vka = vka
        self.sy = sy
        self.z_D = z_D
        self.slope = slope
        self.allow = allow
        self.delr = m.dis.delr.array
        self.delc = m.dis.delc.array
        
        # If not provided all top layer cells will be considered for PA-TBC
        if self.allow is None:
            self.allow = np.full((self.nrow, self.ncol), True)
            
        # Make properties spatial arrays if provided as scalar
        self.hs = self.make_array(self.hs)
        self.A = self.make_array(self.A)
        self.hk = self.make_array(self.hk)
        self.vka = self.make_array(self.vka)
        self.sy = self.make_array(self.sy)
        self.z_D = self.make_array(self.z_D)
        self.slope = self.make_array(self.slope)
        
        # For variable density simulations define density
        self.rho_fresh = 1000
        self.rho_salt = 1025
        
        # Calculations required
        self.vertical_conductance = self.get_vertical_conductance()
    
    
    def get_drn_stress_period_data(self):
        """
        Generate input for the stress_period_data
        dictionary of the DRN package for the given stress period.

        Returns
        -------
        drn_stress_period_data : array (shape depends on MODFLOW version)
            stress period data for each cell holding a DRN. Contains 
            layer, row, column, boundary hydraulic head and conductance for each
            cell. For SEAWAT, bottom elevation of the drain can be included.
            
        """
        
        # Find DRN cells
        is_drn = (self.hghb <= self.z) & self.allow & self.active
        
        # Get row, col and lay index of DRN cells
        row, col = np.where(is_drn)
        lay = self.idx_top_layer[row, col]
        
        # Select respective heads and conductances
        stage = self.z[row, col]
        cond = self.vertical_conductance[row, col]
        
        if self.version == 'mf6':
            drn_stress_period_data = list(zip(zip(lay, row, col), stage, cond))
        elif self.version == 'seawat':
            # Get surface elevation of the boundary condition
            zghb = self.z[row, col]
            drn_stress_period_data = list(zip(lay, row, col, stage, cond, zghb))
            drn_stress_period_data = np.rec.fromrecords(
                drn_stress_period_data,
                dtype = [('k', '<i8'), ('i', '<i8'), ('j', '<i8'),
                         ('stage', '<f4'), ('cond', '<f4'), ('zghb', '<f4')]
            )
        else:
            drn_stress_period_data = list(zip(lay, row, col, stage, cond))
            drn_stress_period_data = np.rec.fromrecords(
                drn_stress_period_data,
                dtype = [('k', '<i8'), ('i', '<i8'), ('j', '<i8'),
                         ('elev', '<f4'), ('cond', '<f4')]
            )
        
        return drn_stress_period_data
    
    
    def get_ghb_stress_period_data(self):
        """
        Generate input for the stress_period_data
        dictionary of the GHB package for the given stress period.

        Returns
        -------
        ghb_stress_period_data : array (shape depends on MODFLOW version)
            stress period data for each cell holding a DRN. Contains 
            layer, row, column, boundary hydraulic head and conductance for each
            cell. For SEAWAT, bottom elevation of the GHB can be included.

        """
        
        # Load regression parameters for empirical correction function
        lm = pd.read_csv(mfpatbc.path_module + 'coef_lm_sim-slopesteady_nwt.csv')
        
        # Saturated aquifer thickness
        D = self.hs - self.z_D
        
        # Calculate phase-averaged intertidal heads
        self.hghb = mfpatbc.pavg_heads(
            self.hs, self.A, self.z, self.T, self.hk, self.sy, D, self.vka, 
            self.slope, lm
        )
        
        # If SEAWAT or MODFLOW-2005 with SWI2 calculate freswater heads
        if (self.version == 'seawat') or (self.version == 'swi2'):
            self.hghb_fresh = self.get_freshwater_heads()
        else:
            self.hghb_fresh = self.hghb
        
        # Find GHB cells
        is_ghb = (self.hghb > self.z) & self.allow & self.active
        
        # Get row, col and lay index of GHB cells
        row, col = np.where(is_ghb)
        lay = self.idx_top_layer[row, col]
        
        # Select respective heads and conductances
        stage = self.hghb_fresh[row, col]
        cond = self.vertical_conductance[row, col]
        
        if self.version == 'mf6':
            ghb_stress_period_data = list(zip(zip(lay, row, col), stage, cond))
        elif self.version == 'seawat':
            # Get surface elevation of the boundary condition
            zghb = self.z[row, col]
            ghb_stress_period_data = list(zip(lay, row, col, stage, cond, zghb))
            ghb_stress_period_data = np.rec.fromrecords(
                ghb_stress_period_data,
                dtype = [('k', '<i8'), ('i', '<i8'), ('j', '<i8'),
                         ('stage', '<f4'), ('cond', '<f4'), ('zghb', '<f4')]
            )
        else:
            ghb_stress_period_data = list(zip(lay, row, col, stage, cond))
            ghb_stress_period_data = np.rec.fromrecords(
                ghb_stress_period_data,
                dtype = [('k', '<i8'), ('i', '<i8'), ('j', '<i8'),
                         ('bhead', '<f4'), ('cond', '<f4')]
            )
            
        return ghb_stress_period_data
    
    
    def get_freshwater_heads(self):
        """
        Calculates fresh water heads according to equation 6 from 
        Post et al. (2007). 
        
        Assumes that GHB reservoir has salt water density, and that the bottom
        of the reservoir is at surface elevation (aquifer-ocean interface).

        Returns
        -------
        hghb_fresh : array of floats (nrow, ncol)
            fresh water heads.
            
        References
        ----------
        Post, V., Kooi, H., & Simmons, C. (2007). 
        Using Hydraulic Head Measurements in Variable-Density Ground Water 
        Flow Analyses. Ground Water, 45(6), 664–671. 
        https://doi.org/10.1111/j.1745-6584.2007.00339.x

        """
        
        hghb_fresh = self.rho_salt / self.rho_fresh * self.hghb -\
            (self.rho_salt - self.rho_fresh) / self.rho_fresh * self.z
            
        return hghb_fresh
    
    
    def get_ssm_stress_period_data(self, conc_ghb = 35):
        """
        Generate input for SSM package of MT3DMS for PA-TBC GHB cells. Defines
        SSM for GHB for all cells with allow = True and active = True, i.e.
        where PA-TBC can be present.
        
        Output is not based on stress period dependent data (high water,
        mean tide level, ...). Thus, is the same for each stress period if
        allow and active do not change and defining first stress period is 
        sufficient.

        Parameters
        ----------
        conc_ghb : float, optional
            Concentration within the GHB reservoir. The default is 35.
            Units as defined in BTN package.

        Returns
        -------
        ssmghb_stress_period_data : Tarray (shape depends on MODFLOW version)
            stress period data for each cell with attributes allow and active True. 
            Contains layer, row, column, concentrationand boundary condition
            identifier for each cell.

        """
        
        # Find all allowed SSM cells for GHB
        is_ssmghb = self.allow & self.active
        
        # Get row, col and lay index of possible GHB cells
        row, col = np.where(is_ssmghb)
        lay = self.idx_top_layer[row, col]
        
        # Get GHB ITYPE for SSM input
        itype = flopy.mt3d.Mt3dSsm.itype_dict()
        
        # Concentration and boundary condition identifier.
        css = np.full(lay.size, conc_ghb)
        itype_ghb = np.full(lay.size, itype['GHB'])
        
        ssmghb_stress_period_data = list(zip(
            lay, row, col, css, itype_ghb
        ))
        
        ssmghb_stress_period_data = np.rec.fromrecords(
            ssmghb_stress_period_data,
            dtype = [('k', '<i8'), ('i', '<i8'), ('j', '<i8'),
                     ('css', '<f4'), ('itype', '<i8')]
        )
        
        return ssmghb_stress_period_data
        
    
    def get_surface_elevation(self):
        """
        Define surface elevations of the aquifer-ocean interface cells.
        
        Row and column combinations with all layers inactive get value of
        upper_bound[nlay-1]. This does not necessarily coincide with the respective
        surface elevation of the aquifer-ocean interface.

        Returns
        -------
        z : array of floats (nrow, ncol)
            Surface elevations.

        """
                
        upper_bound = np.concatenate(
            [self.top.reshape(1, self.nrow, self.ncol), self.botm],
            axis = 0
        )
                
        z = upper_bound[self.lay_top, self.row_top, self.col_top]
        z = z.reshape((self.nrow, self.ncol))
        
        return z


    def get_top_layer(self):
        """
        Find index of layer whose cell top represent surface elevation. Based
        on IBOUND for MODFLOW-2005 and SEAWAT and on IDOMAIN for MODFLOW 6.
        
        All values of IBOUND or IDOMAIN > 0 are considered active.
        The function searches for the first active layer in each row and column
        combination and defines this cell as the aquifer-ocean buondary.
        
        Should idx_layer_top = None in PATBC initialization, it is added as 
        an attribute here. Note, that in this case inactive row and column 
        combinations with no aquifer-ocean interface are set to layer nlay.
        They are not considered when setting up GHB and DRN dictionaries.

        Returns
        -------
        lay_top : array of ints
            layer of cells holding PA-TBC
        row_top : array of ints
            row of cells holding PA-TBC
        col_top : array of ints
            column of cells holding PA-TBC
        active : array of bools (nrow, ncol)
            True, if any layer holds PA-TBC
            

        """
        
        if self.idx_top_layer is None:
            cell_active = self.ibound.copy()
            # Set all positive values (active cell) to 1 for argmax to work correctly.
            cell_active[cell_active > 0] = 1
            
            # Append additional layer to find row column combinations without 
            # active model cells
            cell_active = np.concatenate(
                [cell_active, np.full((1, self.nrow, self.ncol), 1)],
                axis = 0
            )
            self.idx_top_layer = np.argmax(cell_active, axis = 0)
            self.idx_top_layer[self.idx_top_layer == self.nlay] = -1 
        
        # Find first layer
        row_top, col_top = np.where(self.idx_top_layer > -1e30)
        active = self.idx_top_layer > -1
        # Inactive areas are set to layer nlay.
        self.idx_top_layer[~active] = self.nlay - 1
        lay_top = self.idx_top_layer.reshape(-1)
                
        return lay_top, row_top, col_top, active
        
        
    def get_vertical_conductance(self):
        """
        Calculate vertical conductance for the phase-averaged boundary condition.
        Is calculated according to Harbaugh (2005):
            C = vk * A / l
            
            vk: vertical hydraulic conductivity
            
            A: surface of the cell top (representing beach surface)
            
            l: Length along flow path. Here, from surface elevation 
            (elevation of cell's upper boundary) to elevation of
            cell calculation point.

        Returns
        -------
        cond : array of floats (nrow, ncol)
            Vertical conductance of the phase-averaged boundary condition cells.
            
        References
        -----
        Harbaugh, A. W. (2005). MODFLOW-2005, The U.S. Geological Survey Modular 
        Ground-Water Model—the Ground-Water Flow Process (Techniques and Methods 
        No. 6-A16). Reston: U.S. Geological Survey.
        
        """
        
        length_flowpath = 0.5 * (self.z - self.bottom)
        
        DR, DC = np.meshgrid(self.delr, self.delc)
        area = DR * DC
        
        cond = (self.hk * area) / length_flowpath 
        
        return cond
    
    
    def make_array(self, parameter):
        """
        Creates array of size (nrow, ncol) for attribute if attribute is scalar.

        Parameters
        ----------
        parameter : scalar or array
            attribute to check.

        Returns
        -------
        parameter : array (nrow, ncol)
            Array of attribute.

        """
        if not hasattr(parameter, '__len__'):
            parameter = np.full((self.nrow, self.ncol), parameter)
            
        return parameter
    

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


#def define_filepath(m, filepath, )


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
    input of the DRN input file. This is, at least for SEAWAT auxiliar variable 
    IDRNELEV, not possible using the flopy write method for DRN.
    
    See section 'Examples' in the documentation of mfpatbc.PATBC() for the 
    concept on how to use this function.
    
    So far, only supports writing of stress_period_data dictionaries in 
    format required for flopy.modflow.ModflowDrn (MODFLOW-2005 or SEAWAT).
    
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
        print('Start writing DRN input file\n')
        
    # Open file
    f = open(os.path.join(filepath, filename), 'w')
    
    # Line numbering according to Harbaugh (2005)
    # Line 0
    f.write(f'# DRN package, generated by mfpatbc version {mfpatbc.__version__}' +
            f'{newline}')
    
    # Line 2
    f.write(
        f' {mxactb:{precision["cell"]}d} {idrncb:{precision["cell"]}d}' +\
        f'{options_str}{newline}'
    )
    
    # Lines 5 and 6
    write_stress_period_data(f, m.nper, stress_period_data, precision, newline,
                             verbose)
    
    f.close() 
    
    if verbose:
        print('Finished writing DRN input file\n')
 
      
def write_ghb(m, stress_period_data, mxactb = None, ighbcb = 0, options = None,
              filename = None, filepath = None, verbose = False, 
              precision = {'cell': 9, 'data': 15}, newline = '\n'):
    """
    Write General Head Boundary Package (GHB) input file based on stress period 
    data defined e.g. with method PATBC.get_ghb_stress_period_data(). 
    
    This function allows to write auxiliar variables to the stress period data
    input of the GHB input file. This is, at least for SEAWAT auxiliar variable 
    IGHBELEV, not possible using the flopy write method for DRN.
    
    See section 'Examples' in the documentation of mfpatbc.PATBC() for the 
    concept on how to use this function.
    
    So far, only supports writing of stress_period_data dictionaries in 
    format required for flopy.modflow.ModflowGhb (MODFLOW-2005 or SEAWAT).
    
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
        print('Start writing GHB input file\n')
        
    # Open file
    f = open(os.path.join(filepath, filename), 'w')
    
    # Line numbering according to Harbaugh (2005)
    # Line 0
    f.write(f'# GHB package, generated by mfpatbc version {mfpatbc.__version__}' +
            f'{newline}')
    
    # Line 2
    f.write(
        f' {mxactb:{precision["cell"]}d} {ighbcb:{precision["cell"]}d}' +\
        f'{options_str}{newline}'
    )
    
    # Lines 5 and 6
    write_stress_period_data(f, m.nper, stress_period_data, precision, newline,
                             verbose)
    
    f.close()
    
    if verbose:
        print('Finished writing GHB input file\n')
    
    
def write_stress_period_data(f, nper, stress_period_data, precision, newline,
                             verbose):
    """
    Write stress period data in format of stress_period_data dictonaries of
    flopy.modflow stress period based boundary condition packages 
    (MODFLOW-2005 or SEAWAT).
    
    Code follows algorithm for writing transient stress period based boundary
    condition package input files in flopy.utils.utils_list.write_transient.
    
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
    precision : dict
        precision for cell and data information, when writing boundary condition
        input file. 'cell' refers to the layer, row, and column information, 
        'data' refers to stage, conductivity and auxiliary variables. 
    newline : str
        Newline command.
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
            print(f'Write stress period {kper}\n')
        
        # Stress periods before first dictionary entries
        if kper < first_per:
            itmp = 0

        # If kper present in keys write data to file
        elif kper in keys:
            spd_kper = stress_period_data[kper]
            itmp = len(spd_kper)
            
        # If kper not in dictionary keys, reuse previous stress period data
        else:
            itmp = -1
            
        # Line 5
        f.write(f' {itmp:{precision["cell"]}d} {0:{precision["cell"]}d}' +\
                f' # stress period {kper + 1}\n')
        
        # Line 6 if data present
        if kper in keys:
            # Make cell indexing start with 1
            col_names = spd_kper.dtype.names
            for cid in ['k', 'i', 'j']:
                if cid in col_names:
                    spd_kper[cid] += 1
            
            for idx in range(itmp):                
                data = spd_kper[idx]
                
                ndata = len(data)
                
                dtype = ['d', 'd', 'd'] + np.repeat('G', ndata-3).tolist()
                prec = np.repeat(precision['cell'], 3).tolist() + \
                    np.repeat(precision['data'], ndata-3).tolist()
                
                f.write(
                    ''.join(f' {value:{p}{dt}}'\
                        for value, dt, p in zip(data, dtype, prec)) +\
                    f'{newline}'
                )