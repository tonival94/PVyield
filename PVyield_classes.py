# -*- coding: utf-8 -*-
"""
Defines classes used in PVyield
 
@author: Antonio Vallone (-T O N I-)
"""
from pathlib import Path
from dataclasses import dataclass
import warnings
import numpy as np
import pandas as pd
from typing import Tuple, List
from pvlib import tools

# ---------------------------------------------------------------------------
#   CONSTANTS
# ---------------------------------------------------------------------------
_CURRENT_DIR = Path(__file__).parent.resolve() # get current file directory

_PAN_PARAMETERS         = ['Manufacturer', 'Model', 'PNom', 'Height', 'Width', 'NCelS', 'NCelP', 'NDiode', 'Isc', 'Voc', 'Imp', 'Vmp', 'muISC', 'Gamma', 'muGamma', 'RShunt', 'Rp_0', 'Rp_Exp', 'RSerie', 'Absorb', 'BifacialityFactor', 'iam_theta', 'iam_loss', 'VRevDiode', 'SubModuleLayout']
_OND_PARAMETERS         = ['Manufacturer', 'Model', 'PNomConv', 'PMaxOUT', 'VMppMin', 'VMPPMax', 'IMaxDC', 'NbInputs', 'NbMPPT', 'VOutConv', 'TPNom', 'TPMax', 'VNomEff1', 'VNomEff2', 'VNomEff3', 'EfficEuroV1', 'EfficEuroV2', 'EfficEuroV3', 'p_in_VNomEff1', 'p_in_VNomEff2', 'p_in_VNomEff3', 'p_out_VNomEff1', 'p_out_VNomEff2', 'p_out_VNomEff3', 'TPLim1', 'PLim1', 'TPLimAbs','PLimAbs']
_VALID_TMY_SOURCE       = ['pvgis', 'solargis tmy', 'solcast api', 'meteonorm api']
_VALID_MOUNT_TYPE       = ['single-axis tracker', 'single axis tracker', 'tracker', 'fixed tilt', 'fixed', 'fix']
_VALID_MOD_CONFIG       = ['1p','2p','3p','4p','1l','2l','3l','4l','5l','6l']
_VALID_SUB_MODULE       = ['sltwinhalfcells', 'slinlength', 'half cut', 'halfcut', 'full cut', 'fullcut']
_VALID_IAM_LIST         = ['pan', 'user', 'user_defined', 'glassar', 'glass', 'glass ar', 'ar glass', 'glass_ar', 'ar_glass']
_VALID_MV_CABLE_SECT    = [50, 70, 95, 120, 150, 185, 240, 300, 400, 500, 630]  #[mm2]
_VALID_LV_CABLE_SECT    = [50, 70, 95, 120, 150, 185, 240, 300, 400]            #[mm2]
_VALID_CABLE_MATERIAL   = ['alu', 'cu', 'aluminum', 'copper']
_VALID_FIELD_CURTAILMENT_POINT   = [None, 'inverter out', 'field mv cabin']
_VALID_PROJECT_CURTAILMENT_POINT = [None, 'final mv line in', 'hv trafo in', 'hv trafo out', 'hv line out']
_VALID_MV_TREE_EXITS    = ['grid', 'poc', '']
_EMPTY                  = object()

 
# ---------------------------------------------------------------------------
# SIMULATION CONSTANTS
# ---------------------------------------------------------------------------

class Constants:
    """
    Collects every hidden constant used by the PVyield engine.

    These values are not exposed to the user in the launcher because they are either physical constants,
    modelling conventions, or empirical coefficients that should not be changed on a project-by-project basis.
    They are grouped here so that they can be reviewed and tuned in a single place.

    Attributes:
        freq (str)                    [-] timestep of the simulation. Only '1h' is supported, unless a TMY with a finer timestep is imported manually
        sim_year (int)                [-] conventional year used to build the simulation timestamps. Any non leap year works
        horizon_strip_limit (float)   [deg] band of the sky dome where the horizon diffuse component is assumed to come from
        iam_a_r (float)               [-] empirical angular loss coefficient for the Martin-Ruiz IAM model. Standard range is 0.08 to 0.25.
                                          0.105 produces the best match with bank advisors' results. The IAM model is VERY sensitive to this parameter
        I_o_ref (float)               [A] dark (diode reverse saturation) current at reference conditions. 0.018 nA is the standard PVsyst value
        EgRef (float)                 [eV] energy bandgap at reference temperature. 1.121 eV for crystalline silicon. Must be > 0
        irrad_ref (float)             [W/m2] STC irradiance
        temp_ref (float)              [degC] STC cell temperature
        substrings_per_module (int)   [-] number of substrings per module. For half-cut modules it is the number of substrings in HALF a module,
                                          for full-cut modules it is the number of substrings in the WHOLE module
        default_altitude (float)      [m] altitude assumed when the altitude lookup service is not reachable
        default_pressure (float)      [Pa] air pressure assumed when the TMY does not carry pressure data
        surface_to_axis_offset (float)[m] distance between the torque tube and the module surface (for fixed tilt, between support plane and module surface)
        module_efficiency_thermal (float) [-] conventional module efficiency used by the PVsyst thermal model. PVsyst default is 0.1
        max_shading_integration_step (float) [deg] max angular step used to pick the number of rows in the ground view factor integration
        ground_tilt_steps (int)       [-] number of structure tilt angles on which the ground view weights are pre-computed.
                                          The weights are then interpolated at every timestep, the same approach PVsyst uses with
                                          its 7 pre-calculated tracker positions
        ground_npoints (int)          [-] number of ground elements the pitch is discretised into
        ground_n_pvrow (int)          [-] number of points the PV row face is discretised into for the ground view weights
        ground_max_rows (int)         [-] neighbouring rows considered when blocking the lines of sight to the ground
        ground_max_rows_sky (int)     [-] neighbouring rows considered when computing the sky seen from each ground point
        hv_line_loss_per_km (float)   [%] HV line loss per km of line, applied to the HV connection
        resistivity_cu (float)        [ohm*m] copper resistivity, as defined in PVsyst
        resistivity_alu (float)       [ohm*m] aluminium resistivity, as defined in PVsyst
        time_shift_scan (tuple)       [min] (start, stop) of the time shift scan used to align the imported TMY with the simulation timesteps
        solar_constant (float)        [W/m2] solar constant used for the extraterrestrial radiation
    """

    freq                        = '1h'
    sim_year                    = 1990
    horizon_strip_limit         = 6.5
    iam_a_r                     = 0.105
    I_o_ref                     = 0.018E-9
    EgRef                       = 1.121
    irrad_ref                   = 1000
    temp_ref                    = 25.0
    substrings_per_module       = 3
    default_altitude            = 20
    default_pressure            = 101325
    surface_to_axis_offset      = 0.1
    module_efficiency_thermal   = 0.1
    max_shading_integration_step = 5
    ground_tilt_steps           = 15
    ground_npoints              = 100
    ground_n_pvrow              = 15
    ground_max_rows             = 4
    ground_max_rows_sky         = 10
    hv_line_loss_per_km         = 0.1
    resistivity_cu              = 1.68E-8
    resistivity_alu             = 2.82E-8
    time_shift_scan             = (-30, 30)
    solar_constant              = 1361
    numerical_zero              = 1E-6


# ---------------------------------------------------------------------------
# PV YIELD MAIN CLASSES
# ---------------------------------------------------------------------------

@dataclass
class Coordinates:
    """
    Defines Project location via GPS coordinates.

    Attributes:
        latitude (float)   [deg] Latitude, as decimal degrees (e.g. 20.564) - positive is N emisphere, negative is S emisphere [-90, 90]
        longitude (float)  [deg] Longitude, as decimal (e.g. 15.739) - positive is E emisphere, negative is W emisphere. range [-180, 180]
    """

    latitude: float = None
    longitude: float = None
    
    def input_check(self):
        self.latitude = _assign_numeric(self.latitude, 'Project latitude', lo=-90, hi=90)
        self.longitude = _assign_numeric(self.longitude,  'Project longitude', lo=-180, hi=180)
    
@dataclass
class Tmy:
    """
    Defines what TMY source to use for simulation

    Attributes:
        TMY_source (str)              [-] TMY source, either ['PVGIS', 'SolarGIS TMY', 'Solcast API', 'Meteonorm API']. Case sensitive
        api_key (str, optional)       [-] Only for subscription based services like SolarGIS_API, Solcast_API and SolarAnywhere_API. Default to ''.
        TMY_file_name (str, optional) [-] Used only for option 'SolarGIS_tmy'. Save the TMY file in the subfolder named TMY inside the same folder where PVsim is located. Default to ''.
    """

    TMY_source: str
    api_key: str = ''
    TMY_file_name: str = ''
    
    
    def input_check(self):
        self.TMY_source = _assign_str_from_list(self.TMY_source, 'Project TMY_source', _VALID_TMY_SOURCE, standard='pvgis', custom_msg='. TMY will be retrieved automatically from PVGIS web service based on coordinates provided')
        
        # If TMY source need API key, checks that API_key is valid
        if 'API' in self.TMY_source:
            if not self.api_key or self.api_key is None:
                raise ValueError(f'If you wish to import TMY from {self.TMY_source}, you must provide a valid api key')
            else:
                self.api_key = _assign_str(self.api_key, 'api_key')
    
        # If TMY source need TMY file, checks that TMY file path is valid
        if 'TMY' in self.TMY_source:
            if not self.TMY_file_name or self.TMY_file_name is None:
                raise ValueError(f'If you wish to import TMY from {self.TMY_source}, you must provide a valid TMY source')
            else: self.TMY_file_name = _assign_str(self.TMY_file_name, 'TMY_file_name')

@dataclass
class Power:
    """
    Defines DC and AC power

    Attributes:
        p_dc_peak (float)            [MWp] field peak DC power; rounded down to match string count
        p_ac_field_inverter (float)  [MVA] field desired nominal inverter power
    """

    p_dc_peak: float
    p_ac_nom: float
        
    def input_check(self, field_name):
        self.p_dc_peak = _assign_numeric(self.p_dc_peak, field_name+ ' p_dc_peak', lo=0)
        self.p_ac_nom = _assign_numeric(self.p_ac_nom, field_name+ ' p_ac_nom', lo=0)    

@dataclass
class Mms:
    """
    Defines the Module Mounting Structure (MMS) parameters to use in simulation.

    Attributes:
        mount_type (str)               [-] 'Single-axis tracker' or 'Fixed tilt'. NOT case sensitive
        mod_configuration (str)        [-] '1P', '2P', '3P' ... (Portrait) or '1L', '2L', '3L'... (Landscape). NOT case sensitive
        modules_per_string (int)       [-] Number of modules per string in the array.
        pitch (float)                  [m] Row-to-row pitch.
        tilt_angle (float)             [deg] MMS tilt/range angle input. 
        ground_gap (float, optional)   [m] Minimum allowed module-to-ground distance (distance between ground and lower border of PVtable). For tracker, is the gap when MMS is at max tilt.
        axis_azimuth (float, optional) [deg] Compass direction along which the axis of structure lies. Measured clockwise from north (e.g., 90=East, 180=South). If None, inferred from mount type (180 for SAT, 90 for FIX). Default to None.
        backtrack (bool, optional)     [-] SAT parameter - backtracking activation. Default to True.


    About tilt_angle:
    if SAT, tilt_angle variable is the maximum rotation angle, in decimal degrees, of the one-axis tracker from its
    horizontal position (horizontal if axis_tilt = 0). If a float is provided, it represents the maximum rotation angle,
    and the minimum rotation angle is assumed to be the opposite of the maximum angle.
    If a tuple of (min_angle, max_angle) is provided, it represents both the minimum and maximum rotation angles.
    A rotation to 'tilt_angle' is a counter-clockwise rotation about the y-axis of the tracker coordinate system.
    For example, for a tracker with ‘axis_azimuth’ oriented to the south, a rotation to 'tilt_angle' is towards the west,
    and a rotation toward ‘min_angle’ is in the opposite direction, toward the east.
    Hence a max_angle of 180 degrees (equivalent to angle = (-180, 180) allows the tracker to achieve its full rotation capability.

    if FIX, tilt_angle variable is the  surface tilt angles in decimal degrees. The tilt angle is defined as degrees from horizontal
    (e.g. surface facing up = 0, surface facing horizon = 90)
    """
    
    mount_type: str
    mod_configuration: str
    modules_per_string: int
    pitch: float
    tilt_angle: float
    ground_gap: float
    axis_azimuth: float|None = None
    backtrack: bool = False
        
    def input_check(self, field_name):
        self.mount_type = _assign_str_from_list(self.mount_type, field_name+': structure mount type', _VALID_MOUNT_TYPE)
        self.mod_configuration = _assign_str_from_list(self.mod_configuration, field_name+': mod_configuration', _VALID_MOD_CONFIG)
        self.modules_per_string = _assign_numeric(self.modules_per_string, field_name+': modules per string', lo=0, int_only=True)
        self.pitch = _assign_numeric(self.pitch, 'pitch', lo=0)
        
        if 'tracker' in self.mount_type:
            if isinstance (self.tilt_angle, (float, int)):
                self.tilt_angle = (-self.tilt_angle, self.tilt_angle) # transform single int or float value into a tuple to ensure compatibility with all PVlib functions
            self.tilt_angle = tuple(_assign_numeric_array_like(self.tilt_angle, field_name+': tracker tilt angle', exact_size=2))
        else:
            self.tilt_angle = _assign_numeric(self.tilt_angle, field_name+': structure tilt angle')
        
        self.ground_gap = _assign_numeric(self.ground_gap, field_name+': structure ground gap', lo=0, standard=0.5, unit ='m')
 
        if 'tracker' in self.mount_type:
            standard_axis_azimuth = 180
            custom_msg = ', assuming tracker structure is oriented along N-S direction'
            self.backtrack = _assign_bool(self.backtrack, field_name+': backtracking', standard=True)
        else:
            standard_axis_azimuth = 90
            custom_msg = ', assuming fixed structure is oriented along E-W direction'
            self.backtrack = False
            
        
        self.axis_azimuth = _assign_numeric(self.axis_azimuth, field_name+': structure axis azimuth', lo=-360, hi=360, standard=standard_axis_azimuth, unit ='deg', custom_msg=custom_msg )
        
        # if axis_azimuth is negative, a positive value is generated by finding complementary angle. PV lib infact accept only positive angles
        if self.axis_azimuth <0:
            self.axis_azimuth = 360+self.axis_azimuth

@dataclass
class Ground:
    """
    Group all characteristics of ground (albedo &, slope).

    Attributes:
        albedo (float)        [-] Default soil surface albedo. Can be overwritten by imported data (e.g. if using Solcast API)
        slope_tilt (float)    [%] Ground slope in percent
        slope_azimuth (float) [deg] Direction of max ground slope gradient.Is the direction where a ball would roll if left on the slope:
                                    0 = North, 
                                    90 = East,
                                    180 = South
                                    270 = West

    For example, if a site has a slope of 10% in the east direction, the data will be:
    slope tilt = 10
    slope_azimuth = 90
    """
 
    albedo: float
    slope_tilt: float
    slope_azimuth: float
        
    def input_check(self, field_name):
        self.albedo = _assign_numeric(self.albedo, field_name+': albedo', lo=0, hi=1, standard=0.2)
        self.slope_tilt = _assign_numeric(self.slope_tilt, field_name+': slope_tilt', lo=0, standard=0, unit = '%', custom_msg=' (flat surface)')
        self.slope_azimuth = _assign_numeric(self.slope_azimuth, field_name+': slope_azimuth', lo=0, hi=360, standard=0, unit='deg' , custom_msg=' (flat surface)')

@dataclass
class Module:
    """
    Define PV module parameters. 
    
    Can be used in 2 ways:
        - Manually: set PAN_import = False -> User must define all class attributes manually
        - Automatically from PAN: set PAN_import = True and PAN_file_name = 'PAN file path' -> Module class attribute are imported from PAN

    Mandatory attributes:
    -------------------
    PAN_import (bool): True will import class attribtes from PAN file; with False user must define parameters manually.
    
    Attributes to be declared mandatory if PAN_import = True
    -------------------
    PAN_file_name (str, optional):       PAN file name. PAN file must be saved in 'PAN' folder.
        
    Attributes to be declared mandatory if PAN_import = False
    -------------------
    manufacturer (str, optional)        [-] PV module manufacturer.
    model (str, optional)               [-] PV module model.
    Pnom (float, optional)              [Wp] Module DC power at STC.
    height (float, optional)            [m] PV module height.
    width (float, optional)             [m] PV module width.
    cells_in_series (int, optional)     [-] Number of cells connected in series. Assume full-cut cells equivalent (for half-cut, divide actual number by 2).
    cells_in_parallel (int, optional)   [-] Number of cells in parallel. Typically 2 for half-cut and 1 for full-cut.
    n_diodes (int, optional)            [-] Number of bypass diodes, typically 3.
    I_L_ref (float, optional)           [A] short circuit current / Light-generated current (photocurrent) at STC.
    Voc (float, optional)               [V] Open circuit voltage at STC.
    imp_ref (float, optional)           [A] MPP current at STC.
    vmp_ref (float, optional)           [V] MPP voltage at STC.
    alpha_sc (float, optional)          [mAmps/Celsius] Short-circuit current temp coefficient. On PAN is 'muIsc' (e.g. 6.57 for Trina NEG 685).
    gamma_ref (float, optional)         [-] Diode ideality factor. 1.02 is standard PVsyst model.
    mu_gamma (float, optional)          [1/K] Temperature coefficient for the diode ideality factor.
    R_sh_ref (float, optional)          [ohm] Shunt resistance at STC. If None, inferred from other inputs.
    R_sh_0 (float, optional)            [ohm] Shunt resistance at zero irradiance. If None, inferred.
    R_sh_exp (float, optional)          [-] Single diode exponential parameter. Default a 5.5.
    R_s (float, optional)               [ohm] Series resistance at STC. If None, inferred.
    alpha_absorption (float, optional)  [-] Module absorptivity (1 - reflectivity).
    bifaciality (float, optional)       [-] Bifaciality factor.
    theta_ref (list[int], optional)     [degrees] Array of angles at which the IAM is known.
    iam_ref (list[float], optional)     [-] IAM values for each angle in iam_theta.
    V_bypass_diode (float, optional)    [V] Bypass diode voltage drop. Default a 0.7 (per moduli in silicio).
    substring_type (str, optional)      [-] Either 'slHalfCells' (halfcut cell) or 'slInLength' (fullcut cells). NOT case sensitive
    
    Attributes automatically created from PAN if PAN_import = True
    -------------------
    Same list as as above
    """      

    PAN_import: bool
    PAN_file_name: str|None = None
    manufacturer: str|None = None
    model: str|None = None
    Pnom: float|None = None
    height: float|None = None
    width: float|None = None
    cells_in_series: int|None = None
    cells_in_parallel: int|None = None
    n_diodes: float|None = None
    I_L_ref: float|None = None
    Voc: float|None = None
    imp_ref: float|None = None
    vmp_ref: float|None = None
    alpha_sc: float|None = None
    gamma_ref: float|None = None
    mu_gamma: float|None = None
    R_sh_ref: float|None = None
    R_sh_0: float|None = None
    R_sh_exp: float|None = None
    R_s: float|None = None
    alpha_absorption: float|None = None
    bifaciality: float|None = None
    theta_ref: int|None = None
    iam_ref: int|None = None
    V_bypass_diode: float|None = None
    substring_type: str|None = None
    
    def __post_init__(self):
        """
        This post_init function is run automatically right after class initialisation.
        This function import data from PAN file.
        """
        self.PAN_import = _assign_bool(self.PAN_import, 'PAN_import')  # check immediatly that PAN_import is correct value
        
        # If PAN_import = True, import automatically module data from PAN
        if self.PAN_import:
            self.PAN_file_name = _assign_str(self.PAN_file_name, 'PAN_file_name')  # check immediatly that PAN_file_name is a string
            
            PAN_path = _CURRENT_DIR / 'PAN' / self.PAN_file_name  # combines full path of PAN file
            try:
                self.MODULE = _pan_import(file_path = PAN_path,
                                          PAN_parameters = _PAN_PARAMETERS,
                                          pvlib_remap = True,
                                          )
            except Exception as err:
                raise RuntimeError (f'Could not import PV module data from PAN provided. Please make sure PAN file provided is legit and that PAN path provided is correct. Original error: {err}')
        
        # if PAN_import = False, input variables are temporarily stored in a dictionary to mimic what _pan_import function does and to have unified value assignation below
        else:
            self.MODULE = MODULE = {}
            MODULE['manufacturer'] = self.manufacturer
            MODULE['model'] = self.model
            MODULE['Pnom'] = self.Pnom
            MODULE['height'] = self.height
            MODULE['width'] = self.width
            MODULE['cells_in_series'] = self.cells_in_series
            MODULE['cells_in_parallel'] = self.cells_in_parallel
            MODULE['n_diodes'] = self.n_diodes 
            MODULE['I_L_ref'] = self.I_L_ref
            MODULE['Voc'] = self.Voc
            MODULE['imp_ref'] = self.imp_ref
            MODULE['vmp_ref'] = self.vmp_ref
            MODULE['alpha_sc'] = self.alpha_sc
            MODULE['gamma_ref'] = self.gamma_ref
            MODULE['mu_gamma'] = self.mu_gamma
            MODULE['R_sh_ref'] = self.R_sh_ref
            MODULE['R_sh_0'] = self.R_sh_0
            MODULE['R_sh_exp'] = self.R_sh_exp
            MODULE['R_s'] = self.R_s
            MODULE['alpha_absorption'] = self.alpha_absorption
            MODULE['bifaciality'] = self.bifaciality
            MODULE['theta_ref'] = self.theta_ref
            MODULE['iam_ref'] = self.iam_ref
            MODULE['V_bypass_diode'] = self.V_bypass_diode
            MODULE['substring_type'] = self.substring_type
        
    def input_check(self, field_name):
        # check inputs either automatically imported or manually defined
        # we use DICT.get('key') instead of direct access method DICT['key'] because .get methos return a None if key is not found, while direct access retuen an error
        # the None value will then be catched by the assign_type function
        self.manufacturer = _assign_str(self.MODULE.get('manufacturer'), field_name+': PV module manufacturer name', standard= 'User Defined')
        self.model = _assign_str(self.MODULE.get('model'), field_name+': PV model name', standard= 'User Defined')
        self.Pnom = _assign_numeric(self.MODULE.get('Pnom'), field_name+': PV module nominal power PNom', lo = 0)
        self.height = _assign_numeric(self.MODULE.get('height'), field_name+': PV module height', lo = 0)
        self.width = _assign_numeric(self.MODULE.get('width'), field_name+': PV module width', lo = 0)
        self.cells_in_series = _assign_numeric(self.MODULE.get('cells_in_series'), field_name+': PV module cells in series', lo = 0, standard=66, int_only=True, custom_msg=', compatible with recent half-cut modules' )
        self.cells_in_parallel = _assign_numeric(self.MODULE.get('cells_in_parallel'), field_name+': PV module cells in parallel', lo = 0, standard=2, int_only=True, custom_msg=', compatible with recent half-cut modules' )
        self.n_diodes = _assign_numeric(self.MODULE.get('n_diodes'), field_name+': PV module by-pass diode number', lo = 0, standard=3, int_only=True, custom_msg=', which is standard for most of PV modules' )
        self.imp_ref = _assign_numeric(self.MODULE.get('imp_ref'), field_name+': PV module MPPT current imp_ref', lo = 0)
        self.Voc = _assign_numeric(self.MODULE.get('Voc'), field_name+': PV module open circuit voltage Voc', lo = 0)
        self.vmp_ref = _assign_numeric(self.MODULE.get('vmp_ref'), field_name+': PV module MPPT voltage vmp_ref', lo = 0, standard=(self.Pnom/self.imp_ref), unit = 'V', custom_msg=', calculated as Pnom/Imp')    
        self.alpha_sc = _assign_numeric(self.MODULE.get('alpha_sc'), field_name+': PV module short circuit temperature coefficient alpha_sc', lo = 0, standard=6.75, unit = 'mA/Celsius', custom_msg=', which is standard for most of PV modules')
        self.gamma_ref = _assign_numeric(self.MODULE.get('gamma_ref'), field_name+': PV module Gamma (Ideality factor) gamma_ref', lo = 0, standard=1.02, custom_msg=', which is standard for most of PV modules')
        self.mu_gamma = _assign_numeric(self.MODULE.get('mu_gamma'), field_name+': PV module mu_gamma (temperature coefficient for ideality factor) mu_gamma', standard=-0.0002, unit = '1/K', custom_msg=', which is standard for most of PV modules')
        self.I_L_ref = _assign_numeric(self.MODULE.get('I_L_ref'), field_name+': PV module short circuit current I_L_ref', lo = 0, standard=(self.imp_ref + 1), unit = 'A', custom_msg=', derived from other module data') # Standard value formula assume steep V-I curve. Tipically Isc is ~1+Imp
        self.R_sh_ref = _assign_numeric(self.MODULE.get('R_sh_ref'), field_name+': PV module shunt resistance R_sh_ref', lo = 0, standard=(self.vmp_ref/(0.2 * (self.I_L_ref - self.imp_ref))), custom_msg=', derived from other other module data') # Standard value formula is from from PVsyst guide
        self.R_sh_0 = _assign_numeric(self.MODULE.get('R_sh_0'), field_name+': PV module shunt resistance at null irradiance R_sh_0', lo = 0, standard=(self.R_sh_ref*4), unit = 'ohm', custom_msg=', derived from other other module data') # Standard value formula is from from PVsyst guide
        self.R_sh_exp = _assign_numeric(self.MODULE.get('R_sh_exp'), field_name+': PV module exponential factor R_sh_exp', lo = 0, standard=5.5, custom_msg=', which is standard for PV modules single diode model') # Standard value is from from PVsyst guide
        self.R_s = _assign_numeric(self.MODULE.get('R_s'), field_name+': PV module serie resistance exponent R_s', lo = 0, standard=8/(3*self.I_L_ref), unit = 'ohm', custom_msg=', derived from other other module data') # Standard value formula is approximation from datasheet observed
        self.alpha_absorption = _assign_numeric(self.MODULE.get('alpha_absorption'), field_name+': PV module absorption coefficient alpha_absorption', lo = 0, hi=1, standard=0.9, custom_msg=', which is standard for most of PV modules') # Standard value is from from PVsyst guide
        self.bifaciality = _assign_numeric(self.MODULE.get('bifaciality'), field_name+': PV module bifaciality factor', lo = 0, hi=1, standard=0.8, custom_msg=', which is standard for most of PV modules') # Standard value is from from PVsyst guide
        self.theta_ref = self.MODULE.get('theta_ref')
        self.iam_ref = self.MODULE.get('iam_ref')
        self.V_bypass_diode = _assign_numeric(self.MODULE.get('V_bypass_diode'), field_name+': PV module by-pass diode voltage', standard = 0.7, custom_msg=', which is standard for most of PV modules')
        self.V_bypass_diode = abs(self.V_bypass_diode)
        self.substring_type = _assign_str_from_list(self.MODULE.get('substring_type'), field_name+': PV module substring layout', _VALID_SUB_MODULE, standard = 'half cut', custom_msg=' (half-cut cell module)')
        
        del self.MODULE           
        """
        About R_s calculation
        Rserie is the slope of the IV curve at Voc. In an ideal diode is 0
        It can be overestimated as (Voc-Vmp)/(Imp) (slope between Voc and Max power point)
        
        To match the slope, I found out that this formula works better: (Voc-Vmp)/(3xImp)
        If Voc is not provided, then we assume that Voc-Vmp = 8 V 
        """
        
        _compare_numeric(self.imp_ref, self.I_L_ref, field_name+' - PV module max power point current', field_name+' - PV module short circuit current', custom_msg='. Please check PAN file or PV module inputs of '+field_name)
        _compare_numeric(self.vmp_ref, self.Voc, field_name+' - PV module max power point voltage', field_name+' - PV module open circuit', custom_msg='. Please check PAN file or PV module inputs of '+field_name)
        

@dataclass
class Inverter:
    """
    Define Inverter module parameters. 
    
    Can be used in 2 ways:
        - Manually: set OND_import = False -> User must define all class attributes manually
        - Automatically from OND: set OND = True and OND_file_name = 'OND file path' -> Module class attribute are imported from OND

    Mandatory attributes:
    -------------------
    OND_import (bool): [-] True will import class attribtes from OND file; with False user must define parameters manually.
    
    Attributes to be declared if OND_import = True
    -------------------
    OND_file_name (str, optional):    [-] OND file name. OND file must be saved in 'OND' folder.
        
    Attributes to be declared if OND_import = False
    -------------------
    manufacturer (str, optional)     [-] PV module manufacturer.
    model (str, optional)            [-] PV module model.
    Pnom_ac (float, optional)        [kW] nominal AC output power for inverter
    P_max_ac (float, optional)       [kW] max AC output power for inverter
    v_min_mppt  (float, optional)    [V] min MPPT voltage to allow MPPT to work 
    v_max_mppt (float, optional)     [V] max MPPT voltage to allow MPPT to work
    i_max_mppt (float, optional)     [A] max current per MPPT input
    n_string_inputs (int, optional)  [-] number of maximum string per inverter
    n_mppt (int, optional)           [-] number of inverter MPPT
    v_out (float, optional)          [V] inverter output voltage
    inv_euro_eff (float, optional)   [%] average inverter efficiency (EURO efficiency from datasheet)
    T_p_nom (float, optional)        [°C] outdoor temperature at which the inverter produces the nominal power - default: 40
    T_p_max (float, optional)        [°C] outdoor temperature at which the inverter produces the max power - default: 30
    
    Attributes automatically created from OND if OND_import = True
    -------------------    
    Same as above, excluding inv_euro_eff, which is replaced by more effective variables below:
        
    v_eff1 (float)               [V] inverter output voltage at which the efficiency is known for p_in_v_eff1 and p_out_v_eff1
    v_eff2 (float)               [V] inverter output voltage at which the efficiency is known for p_in_v_eff2 and p_out_v_eff2
    v_eff3 (float)               [V] inverter output voltage at which the efficiency is known for p_in_v_eff3 and p_out_v_eff3
    p_in_v_eff1 (list[float])   [kWp] inverter input DC power producing an putput power of p_out_v_eff1 at v_eff1. Expressed as a list, each item of the list correspond to respective item in p_out_v_eff1
    p_in_v_eff2 (list[float])   [kWp] inverter input DC power producing an putput power of p_out_v_eff2 at v_eff2. Expressed as a list, each item of the list correspond to respective item in p_out_v_eff2
    p_in_v_eff3 (list[float])   [kWp] inverter input DC power producing an putput power of p_out_v_eff3 at v_eff3. Expressed as a list, each item of the list correspond to respective item in p_out_v_eff3
    p_out_v_eff1 (list[float])  [kWac] inverter output AC power resulting from input power of  p_in_v_eff1 at v_eff1. Expressed as a list, each item of the list correspond to respective item in p_in_v_eff1
    p_out_v_eff2 (list[float])  [kWac] inverter output AC power resulting from input power of  p_in_v_eff2 at v_eff2. Expressed as a list, each item of the list correspond to respective item in p_in_v_eff2
    p_out_v_eff3 (list[float]   [kWac] inverter output AC power resulting from input power of  p_in_v_eff3 at v_eff3. Expressed as a list, each item of the list correspond to respective item in p_in_v_eff3
    T_lim_1 (float)             [°C] outdoor temperature that derates  inverter power to PLim1
    T_lim_2 (float)             [°C] outdoor temperature that derates  inverter power to PLimAbs
    p_out_T_lim_1 (float)       [kWac] inverter max derated output power output when outdoor temperature is equal to T_lim_1
    p_out_T_lim_2 (float)       [kWac] inverter max derated output power output when outdoor temperature is equal to T_lim_2
    """    
     
    OND_import: bool
    OND_file_name: str|None = None
    manufacturer: str|None = None
    model: str|None = None
    Pnom_ac: float|None = None
    P_max_ac: float|None = None
    v_min_mppt: float|None = None
    v_max_mppt: float|None = None
    i_max_mppt: float|None = None
    n_string_inputs: float|None = None
    n_mppt: float|None = None
    v_out: float|None = None
    inv_euro_eff: float|None = None
    T_p_nom: float|None = None
    T_p_max: float|None = None
        
    def __post_init__(self):
        """
        This post_init function is run automatically right after class initialisation.
        This function import data from OND file.
        """
        self.OND_import = _assign_bool(self.OND_import, 'OND_import')  # check immediatly that OND_import is correct value
        
        # If OND_import = True, import automatically inverter data from OND
        if self.OND_import:
            self.OND_file_name = _assign_str(self.OND_file_name, 'OND_file_name')  # check immediatly that OND_file_name is a string
            
            OND_path = _CURRENT_DIR / 'OND' / self.OND_file_name  # combines full path of PAN file
            try:
                self.INVERTER = _ond_import(file_path = OND_path,
                                            OND_parameters = _OND_PARAMETERS,
                                            pvlib_remap = True
                                            )
                self.INVERTER['i_max_mppt'] = self.INVERTER.get('i_dc_max')/self.INVERTER.get('n_mppt')
            except Exception as err:
                raise RuntimeError (f'Could not import PV module data from OND provided. Please make sure OND file provided is legit and that OND path provided is correct. Original error: {err}')
        
        
        # if OND_import = False, input variables are temporarily stored in a dictionary to mimic what _ond_import function does and to have unified value assignation below
        else:
            self.INVERTER = INVERTER = {}
            INVERTER['manufacturer'] = self.manufacturer
            INVERTER['model'] = self.model
            INVERTER['Pnom_ac'] = self.Pnom_ac
            INVERTER['P_max_ac'] = self.P_max_ac
            INVERTER['v_min_mppt'] = self.v_min_mppt
            INVERTER['v_max_mppt'] = self.v_max_mppt
            INVERTER['i_max_mppt'] = self.i_max_mppt
            INVERTER['n_string_inputs'] = self.n_string_inputs
            INVERTER['n_mppt'] = self.n_mppt
            INVERTER['v_out'] = self.v_out
            INVERTER['inv_euro_eff'] = self.inv_euro_eff
            INVERTER['T_p_nom'] = self.T_p_nom
            INVERTER['T_p_max'] = self.T_p_max
            
    def input_check(self, field_name):
        # check inputs either automatically imported or manually defined
        # we use DICT.get('key') instead of direct access method DICT['key'] because .get methos return a None if key is not found, while direct access retuen an error
        # the None value will then be catched by the assign_type function
        self.manufacturer = _assign_str(self.INVERTER.get('manufacturer'), 'Inverter manufactuer name', standard= 'User Defined')
        self.model = _assign_str(self.INVERTER.get('model'), 'Inverter model name', standard= 'User Defined')
        self.Pnom_ac = _assign_numeric(self.INVERTER.get('Pnom_ac'),'Inverter nominal AC power', lo=0)
        self.P_max_ac = _assign_numeric(self.INVERTER.get('P_max_ac'),'Inverter max AC power', lo=0)
        self.v_min_mppt = _assign_numeric(self.INVERTER.get('v_min_mppt'),'Inverter MPPT minimum voltage', lo=0)
        self.v_max_mppt = _assign_numeric(self.INVERTER.get('v_max_mppt'),'Inverter MPPT maximum voltage', lo=0)
        self.n_string_inputs = _assign_numeric(self.INVERTER.get('n_string_inputs'), 'Number of string input of inverter', lo=0, int_only=True)
        self.n_mppt = _assign_numeric(self.INVERTER.get('n_mppt'), 'Number of MPPT of inverter', lo=0, int_only=True)
        self.i_max_mppt = _assign_numeric(self.INVERTER.get('i_max_mppt'),'inverter max DC current per MPPT', lo=0)
        self.v_out = _assign_numeric(self.INVERTER.get('v_out'),'Inverter AC output voltage', lo=0, standard=800, unit='V', custom_msg=', compatible with modern string, large scale inverter')
        # The inverter delivers its MAXIMUM power up to T_p_max and its NOMINAL power up to T_p_nom, then it derates.
        self.T_p_nom = _assign_numeric(self.INVERTER.get('T_p_nom'),'Inverter temperature for nominal power', lo=0, standard=40, unit='°C', custom_msg=', compatible with modern string, large scale inverter')
        self.T_p_max  = _assign_numeric(self.INVERTER.get('T_p_max'),'Inverter temperature for maximum power', lo=0, standard=30, unit='°C', custom_msg=', compatible with modern string, large scale inverter')
        
        if not self.OND_import:
            self.inv_euro_eff = _assign_numeric(self.INVERTER.get('inv_euro_eff'),'Inverter average european efficiency', lo=0, hi=100, standard=98.8, unit='%', custom_msg=', compatible with modern string, large scale inverter')
        
        else:
            self.v_eff1 = _assign_numeric(self.INVERTER.get('v_eff1'), 'Inverter efficiency voltage 1 (from OND)', lo = 0)
            self.v_eff2 = _assign_numeric(self.INVERTER.get('v_eff2'), 'Inverter efficiency voltage 2 (from OND)', lo = 0)
            self.v_eff3 = _assign_numeric(self.INVERTER.get('v_eff3'), 'Inverter efficiency voltage 3 (from OND)', lo = 0)
            self.p_in_v_eff1 = _assign_numeric_array_like(self.INVERTER.get('p_in_v_eff1'), 'Inverter power vs efficiency - voltage 1 (from OND)', lo = 0)
            self.p_in_v_eff2 = _assign_numeric_array_like(self.INVERTER.get('p_in_v_eff2'), 'Inverter power vs efficiency - voltage 2 (from OND)', lo = 0)
            self.p_in_v_eff3 = _assign_numeric_array_like(self.INVERTER.get('p_in_v_eff3'), 'Inverter power vs efficiency - voltage 3 (from OND)', lo = 0)
            self.p_out_v_eff1 = _assign_numeric_array_like(self.INVERTER.get('p_out_v_eff1'), 'Inverter power vs efficiency - voltage 1 (from OND)', lo = 0)
            self.p_out_v_eff2 = _assign_numeric_array_like(self.INVERTER.get('p_out_v_eff2'), 'Inverter power vs efficiency - voltage 2 (from OND)', lo = 0)
            self.p_out_v_eff3 = _assign_numeric_array_like(self.INVERTER.get('p_out_v_eff3'), 'Inverter power vs efficiency - voltage 3 (from OND)', lo = 0)
            self.T_lim_1 = _assign_numeric(self.INVERTER.get('T_lim_1'), 'Inverter derating temperature 1 (from OND)', lo = 0)
            self.T_lim_2 = _assign_numeric(self.INVERTER.get('T_lim_2'), 'Inverter derating temperature 2 (from OND)', lo = 0)
            self.p_out_T_lim_1 = _assign_numeric(self.INVERTER.get('p_out_T_lim_1'), 'Inverter derated AC power 1 (from OND)', lo = 0)
            self.p_out_T_lim_2 = _assign_numeric(self.INVERTER.get('p_out_T_lim_2'), 'Inverter derated AC power 2 (from OND)', lo = 0)
            del self.INVERTER  
            
            _compare_numeric(self.T_lim_1, self.T_lim_2, 'Inverter power derating temperature 1', 'Inverter power derating temperature 2', custom_msg='. Check OND file')
            _compare_numeric(self.p_out_T_lim_2, self.p_out_T_lim_1, 'Inverter power derated 2', 'Inverter power derated 1', custom_msg='. Check OND file')
     
            
        # compare values to ensure smaller/bigger values are respected
        _compare_numeric(self.Pnom_ac, self.P_max_ac, field_name+' - Inverter nominal AC power', field_name+' - inverter max AC power', custom_msg='. Check OND file or inverter inputs')
        _compare_numeric(self.v_min_mppt, self.v_max_mppt, field_name+' - Inverter MPPT minimum voltage', field_name+' - inverter MPPT maximum voltage', custom_msg='. Check OND file or inverter inputs')   
        _compare_numeric(self.T_p_max, self.T_p_nom, field_name+' - Inverter temperature at maximum power', field_name+' - Inverter temperature at nominal power', custom_msg='. Check OND file or inverter inputs')

@dataclass
class MV_Trafo:
    """
    Defines the LV/MV transformer parameters to use in simulation.

    Attributes:
        v_out (float)                    [kV] LV/MV trafo output voltage (same as connection voltage if no HV connection is present). Tipically 15, 20 or 30kV
        p_max (float, optional)          [MVA] LV/MV transformer max rated power - the simulation will adapt power to field need
        iron_loss (float, optional)      [%] iron loss in the LV/MV transformer - default is 0.1%
        copper_loss (float, optional)    [%] copper loss in the LV/MV transformer - default is 1.0%
   """ 
    v_out: float
    p_max: float|None = None
    iron_loss: float|None = None
    copper_loss: float|None = None
   
    def input_check(self, field_name):
        self.v_out = _assign_numeric(self.v_out, field_name+': LV/MV transformer output voltage', lo=0)
        self.p_max = _assign_numeric(self.p_max, field_name+': LV/MV transformer max rated power', lo=0, standard= 6.4, unit= 'MVA')
        self.iron_loss = _assign_numeric(self.iron_loss, field_name+': LV/MV transformer iron loss', lo=0, hi= 100, standard=0.1, unit='%', custom_msg=', which is standard assumption')
        self.copper_loss = _assign_numeric(self.copper_loss, field_name+': LV/MV transformer copper loss', lo=0, hi= 100, standard=1.0, unit='%', custom_msg=', which is standard assumption')
        
        # compare values to ensure smaller/bigger values are respected
        _compare_numeric(self.iron_loss, self.copper_loss, field_name+'- LV/MV trafo iron loss', field_name+'- LV/MV trafo copper loss', custom_msg='. Check MV trafo inputs input')

@dataclass
class Field_Curtailment:
    """
    Defines field curtailment 

    Attributes:
        grid_curtail_point  (optional,str)   [-] where apply grid curtail. Available options are: (NOT case sensitive)
                                                     # None              : No grid curtailment                                        
                                                     # 'Inverter out'    : the grid curtailment is applied at inverter output. AC output power of the inverter is limited to Power.p_ac_nom for each PV_field. Useful when inverter is curtailed
                                                     # 'Field MV cabin'  : the grid curtailment is applied at MV cabin of the field. Useful when PV Field CCI/PPC curtailment is used or PoC is at PV Field collection cabin
    grid_curtail_power (optional, float)    [MVA] grid AC power curtail for point specified
    """
    grid_curtail_point: str|None = None
    grid_curtail_power: float|None = None
    
    def input_check(self, field_name):
            # the value is normalised (lowercase, single spaces) and remapped through the alias table, so that user input is NOT case sensitive
            self.grid_curtail_point = _assign_str_from_list(self.grid_curtail_point, field_name+': grid_curtail_point', _VALID_FIELD_CURTAILMENT_POINT, standard = None, custom_msg=' hence no curtailment applied')
            if self.grid_curtail_point is not None:
                self.grid_curtail_power = _assign_numeric(self.grid_curtail_power, field_name+': grid_curtail_power', lo=0)

@dataclass
class Loss:
    """
    Defines simulation losses

    Attributes:
        iam_model (str, optional)               [-] 'PAN' or 'GlassAR' (NOT case sensitive)
                                                      'PAN' will use variables iam_theta and iam_loss in Module object
                                                      'GlassAR' will use Martin-ruiz model adjusted for Anti-reflective coated glass
        soiling_loss (float, optional)          [%] soiling loss - default is 2%
        rear_shade_factor (float, optional)     [%] rear shade factor (fraction of light not hitting the back due to bus and bars). 5% is standard PVsyst assumptions
        thermal_loss_factor (float, optional)   [W/m2/K] thermal loss factor for PV module. Default Value is 29, tipical value for most ground mounted applications
        module_quality (float, optional)        [%] module quality loss (a negative value is a gain, a positive value is a loss) - default is -0.3%
        lid_letid_loss (float, optional)        [%] LID + LeTID loss - default is 0.5%
        first_year_degrad (float, optional)     [%] First year degradation - default is 0.3%
        mismatch_loss_front (float, optional)   [%] array mismatch loss for front cells - default is 0.5%
        mismatch_loss_back (float, optional)    [%] array mismatch loss for back cells - default is 0.2%
        dc_ohmic_loss_STC (float, optional)     [%] DC cabling ohmic loss in STC condition - default is 1.5%
        auxiliaries_loss (float, optional)      [%] Auxiliary loss. If left empty, '', or None, it will be inferred from MMS type
        pv_availability (float, optional)       [%] PV availability - if None, will infer the value from MMS: 99% for SAT, 99.5 for FT
        ac_lv_cable_len (float, optional)       [m] average lenght of the AC LV cable from inverter to trafo
        ac_lv_cable_sect (float, optional)      [mm2] average section of the AC LV cable
        ac_lv_cable_material (float, optional)  [-] AC LV cable material. Either 'Cu' (copper) or 'Alu' (aluminum). NOT case sensitive
        ac_mv_cable_len (float, optional)       [m] average lenght of the AC MV cable from trafo to collection Cabin
        ac_mv_cable_sect (float, optional)      [mm2] average section of the AC MV cable from trafo to collection Cabin
        ac_mv_cable_material (float, optional)  [-] AC MV cable material. Either 'Cu' (copper) or 'Alu' (aluminum). NOT case sensitive
    """
    
    iam_model: str = ''
    soiling_loss: float|None = None
    rear_shade_factor: float|None = None
    thermal_loss_factor: float|None = None
    module_quality: float|None = None
    lid_letid_loss: float|None = None
    first_year_degrad: float|None = None
    mismatch_loss_front: float|None = None
    mismatch_loss_back: float|None = None
    dc_ohmic_loss_STC: float|None = None
    auxiliaries_loss: float|None = None
    pv_availability: float|None = None
    ac_lv_cable_len: float|None = None
    ac_lv_cable_sect: float|None = None
    ac_lv_cable_material: str = ''        # here we assign default value so that str.lower will not crash by accepting None as value
    ac_mv_cable_len: float|None = None
    ac_mv_cable_sect: float|None = None
    ac_mv_cable_material: str = ''        # here we assign default value so that str.lower will not crash by accepting None as value
    
        
    def input_check (self, field_name):
        self.iam_model = _assign_str_from_list(self.iam_model, field_name+': iam_model', _VALID_IAM_LIST, standard='glass_ar', custom_msg=', (glass with Anti Reflective coating), which is most conservative')
        self.soiling_loss = _assign_numeric(self.soiling_loss, field_name+': soiling_loss', lo=0, hi=100, standard=2, unit='%', custom_msg=', compatible with most soiling')
        self.rear_shade_factor = _assign_numeric(self.rear_shade_factor, field_name+': rear_shade_factor', lo=0, hi=100, standard=5, unit='%', custom_msg=', compatible with most PV structure')
        self.thermal_loss_factor = _assign_numeric(self.thermal_loss_factor, field_name+': thermal_loss_factor', lo=0, standard=29, unit='W/m2/K', custom_msg=', which is standard assumption for PV modules on ground structure (either fix or tracker)')
        self.module_quality = _assign_numeric(self.module_quality, field_name+': module_quality', lo=-100, hi= 100, standard=-0.3, unit='%', custom_msg=', which is standard assumption')
        self.lid_letid_loss = _assign_numeric(self.lid_letid_loss, field_name+': lid_letid_loss', lo=0, hi= 100, standard=0.5, unit='%', custom_msg=', which is standard assumption for N-type PV modules')
        self.first_year_degrad = _assign_numeric(self.first_year_degrad, field_name+': first_year_degrad', lo=0, hi= 100, standard=0.3, unit='%', custom_msg=', which is standard assumption for N-type PV modules')
        self.mismatch_loss_front = _assign_numeric(self.mismatch_loss_front, field_name+': mismatch_loss_front', lo=0, hi= 100, standard=0.5, unit='%', custom_msg=', which is standard assumption')
        self.mismatch_loss_back = _assign_numeric(self.mismatch_loss_back, field_name+': mismatch_loss_back', lo=0, hi= 100, standard=0.2, unit='%', custom_msg=', which is standard assumption')
        self.dc_ohmic_loss_STC = _assign_numeric(self.dc_ohmic_loss_STC, field_name+': dc_ohmic_loss_STC', lo=0, hi= 100, standard=1.5, unit='%', custom_msg=', which is standard assumption')
        # auxiliaries_loss is assigned in PV_field.input_check
        # pv_availability is assigned in PV_field.input_check
        self.ac_lv_cable_len = _assign_numeric(self.ac_lv_cable_len, field_name+': ac_lv_cable_len', lo=0, standard=500, unit='m', custom_msg=', which is compatible with string inverters')
        self.ac_lv_cable_sect = _assign_numeric_from_list(self.ac_lv_cable_sect, field_name+': ac_lv_cable_sect',_VALID_LV_CABLE_SECT, standard=400, unit='mm2', custom_msg=', which is compatible with 300kVA string inverters')
        self.ac_lv_cable_material = _assign_str_from_list(self.ac_lv_cable_material, field_name+': ac_lv_cable_material', _VALID_CABLE_MATERIAL, standard='copper', custom_msg=' cables')
        self.ac_mv_cable_len = _assign_numeric(self.ac_mv_cable_len, field_name+': ac_mv_cable_len', lo=0, standard=1000, unit='m', custom_msg=', compatible with projects up to 15/20MW')
        self.ac_mv_cable_sect = _assign_numeric_from_list(self.ac_mv_cable_sect, field_name+': ac_mv_cable_sect', _VALID_MV_CABLE_SECT, standard=630, unit='mm2', custom_msg=', which is max commercial cable size')
        self.ac_mv_cable_material = _assign_str_from_list(self.ac_mv_cable_material, field_name+': ac_mv_cable_material', _VALID_CABLE_MATERIAL, standard='aluminum', custom_msg=' cables')      

@dataclass
class MV_Connection:
    """
    Defines MV connection parameters
    
    Attributes:
        line_len (optional, float)    [km] lenght of the MV connection cable, from collection cabin to HV trafo (if present) or PoC. If None, will assume no MV connection
        material (optional, str)      [str] material of MV connection, either 'Cu' (copper) or 'Alu' (aluminum). Not case sensitive
        destination_field (str)       [str] name of the destination field where the connection is going. Must match an existing PV_Field.name or be 'PoC'/'grid'
        auto_sizing (optional, bool)  [bool] if true, number of cable and section will be calculated automatically
        line_nr (optional,float)      [-] number of MV connection cable, from  collection cabin to HV trafo (if present) or PoC.
        line_sect (optional,float)    [mm2] section of the MV connection cable, from collection cabin to HV trafo (if present) or PoC
        """
        
    line_len: float|None = None
    material: str|None = None
    destination_field: str = 'grid'
    auto_sizing: bool = True
    line_nr: int|None = None
    line_sect: float|None = None
                
    
    def input_check(self, field_name):
        self.line_len = _assign_numeric(self.line_len, field_name+': line lenght', lo=0, standard=None, unit='km', custom_msg=' (no MV connection line)')
        
        # if #line_len is not None, assign all other class attributes
        if self.line_len is not None:
            self.material = _assign_str_from_list(self.material, field_name+': line material', _VALID_CABLE_MATERIAL, standard='aluminum', custom_msg=' cables')
            self.auto_sizing = _assign_bool(self.auto_sizing, field_name+': cable sizing', standard = True, custom_msg= ', thus MV cable size and number will be automatically calculated')
            self.destination_field = _assign_str(self.destination_field, field_name+': destination field')

            if not self.auto_sizing:
                self.line_nr = _assign_numeric(self.line_nr, field_name+': line_nr', lo=0, unit='', int_only = True)
                self.line_sect = _assign_numeric_from_list(self.line_sect, field_name+': line_sect', _VALID_MV_CABLE_SECT, standard=630, unit='mm2', custom_msg=', which is max commercial size')      

        """
            else:
                 #!!! inserire qui autosizing!!!
        """

@dataclass
class PV_Field:
    """
    Groups all parameters of a PV field
    
    Attributes:
        name (string)                           [-] name of PV field
        power (Power)                           object defining DC and AC power for PV field
        ground (Ground)                         object defining ground orientation of PV field
        mms (Mms)                               object defining MMS used in PV field          
        module (Module)                         object defining PV modules used in PV field
        inverter (Inverter)                     object defining inverter used in PV field
        mv_trafo (MV_Trafo)                     object defining MV trafo inside PV field
        field_curtailment (Field_Curtailment)   object defining if the re if curtailment at PV field level 
        loss (Loss)                             object defining PV field loss
        mv_connection (MV_Connection)           object defining MV connection outside PV field
        incoming_fields (list[str])             incoming MV_Connection in this PV_Field
        destination_field (str)                 pv field where this field is connetec to via MV_Connection   
    """
    
    name: str
    power: Power
    mms: Mms
    ground: Ground
    module: Module
    inverter: Inverter
    mv_trafo: MV_Trafo
    mv_connection: MV_Connection
    field_curtailment : Field_Curtailment
    loss: Loss
    
    def __post_init__(self):
        # Creates list to memorise all PV_Fields connected to this one via MV_connection.
        # The list will be populated by PV_Project class method _build_mv_tree(self)
        self.incoming_fields = []
        #retrieve from MV_connection the destination field
        self.destination_field = self.mv_connection.destination_field
    
    def input_check(self):
        # 0 - refresh the destination field from the MV_Connection object.
        # __post_init__ takes a snapshot of it at construction time, so a MV_Connection edited afterwards
        # (e.g. re-routing a field while building the project) would otherwise leave a stale destination
        # behind and the MV tree would be built on the old topology.
        self.destination_field = self.mv_connection.destination_field

        # 1 - check minimum pitch
        if 'tracker' in  self.mms.mount_type:
            angle = 0
        else:
            angle = self.mms.tilt_angle
        
        if self.mms.mod_configuration [1] =='p':
            min_pitch = self.module.height * int(self.mms.mod_configuration[0]) * tools.cosd(angle)
        elif self.mms.mod_configuration [1] =='l':
            min_pitch = self.module.width * int(self.mms.mod_configuration[0]) * tools.cosd(angle)
        
        if self.mms.pitch < min_pitch:
            raise ValueError (f'PV field "{self.name}": the pitch you selected is incompatible with the configuration and module size. The pitch should be at least {min_pitch} m')
       
        
        # 2 - check IAM loss compatibility
        if self.loss.iam_model in ['pan', 'user', 'user defined', 'user_defined']:
            if not isinstance (self.module.theta_ref, np.ndarray) or not isinstance (self.module.iam_ref, np.ndarray):
                warnings.warn(f'PV field "{self.name}": Invalid or missing IAM loss values. Loss values should be a list. Simulation will calculate IAM loss assuming Anti-Reflective glass')
                self.loss.iam_model = 'glass_ar'
            elif len(self.module.theta_ref) != len(self.module.iam_ref):
                warnings.warn(f'PV field "{self.name}": Invalid or missing IAM loss values. Loss values should be a list. Simulation will calculate IAM loss assuming Anti-Reflective glass')
                self.loss.iam_model = 'glass_ar'
        if self.loss.iam_model == 'pan' and not self.module.PAN_import:
            warnings.warn(f'PV field "{self.name}": Incompatible IAM loss method. "PAN" method was selected, but no PAN file used in simulation. Simulation will calculate IAM loss assuming Anti-Reflective glass')
            self.loss.iam_model = 'glass_ar'
        
        # 3 - check loss based on MMS
        if 'tracker' in self.mms.mount_type:
            standard_aux_loss = 0.4
            standard_pv_availability = 99
        else: 
            standard_aux_loss = 0.2
            standard_pv_availability = 99.5
            
        self.loss.auxiliaries_loss = _assign_numeric(self.loss.auxiliaries_loss, self.name+': auxiliaries_loss', lo=0, hi=100, standard = standard_aux_loss, unit='%', custom_msg=f', compatible with {self.mms.mount_type} mount type')     
        self.loss.pv_availability = _assign_numeric(self.loss.pv_availability, self.name+': pv_availability', lo=0, hi=100, standard = standard_pv_availability, unit='%', custom_msg=f', compatible with {self.mms.mount_type} mount type')  
        
        
        # 4 - Check module and inverter compatibility
        min_string_len = np.ceil(self.inverter.v_min_mppt/self.module.vmp_ref)
        max_string_len = np.floor(self.inverter.v_max_mppt/self.module.Voc)
        max_dc_ac_ratio = (self.module.Pnom * self.mms.modules_per_string * self.inverter.n_string_inputs) / (self.inverter.Pnom_ac * 1000)
        
        if self.mms.modules_per_string < min_string_len or self.mms.modules_per_string > max_string_len:
            raise RuntimeError(f'PV field "{self.name}": PV modules per string not compatible with inverter. Modules per string must be between {min_string_len} and {max_string_len}')
        if self.module.I_L_ref > self.inverter.i_max_mppt:
            raise RuntimeError(f'PV field "{self.name}": the PV module power is too high for the inverter. Consider using a PV module with lower rated power or an inverter with higher MPPT current limit')
        
        
        # 5 - compute precise DC and AC power
        self.power.p_dc_peak = round(self.power.p_dc_peak*1000*1000/self.module.Pnom/self.mms.modules_per_string)*self.module.Pnom*self.mms.modules_per_string/1000/1000
        warnings.warn(f'PV field "{self.name}": the field DC power was rounded to {self.power.p_dc_peak} MWp to match an integer number of strings')
        
        n_inv_field = np.ceil(self.power.p_ac_nom*1000/self.inverter.Pnom_ac)
        self.power.p_ac_nom = n_inv_field * self.inverter.Pnom_ac /1000
        warnings.warn(f'PV field "{self.name}": the field AC power was rounded to {self.power.p_ac_nom} MVA to match an integer number of inverter')
        
        # 6 - Check DC/AC ratio compatibility
        if self.power.p_dc_peak/self.power.p_ac_nom > max_dc_ac_ratio:
            raise RuntimeError(f'PV field "{self.name}": the desired DC/AC ratio is not theoretically achievable with the chosen module, string length and inverter. Reduce the DC/AC ratio or increase module per string') 

        # 7 - Check field curtailment
        # It is now compared against max_field_ac_power, i.e. the sum of the max AC power of every inverter in the field, in MVA.
        if self.field_curtailment.grid_curtail_point is not None:
            max_field_ac_power = np.ceil(self.power.p_ac_nom*1000/self.inverter.Pnom_ac)*self.inverter.P_max_ac/1000
            if self.field_curtailment.grid_curtail_power > max_field_ac_power:
                warnings.warn(f'PV field "{self.name}": AC power curtailment is ({self.field_curtailment.grid_curtail_power}) MVA but should be lower than the sum of maximum AC power of all inverters ({max_field_ac_power} MVA)')

@dataclass
class HV_Trafo:
    """
    Defines the MV/HV transformer parameters to use in simulation.

    Attributes:
        p_nom (float, optional)          [MVA] MV/HV transformer rated power
        iron_loss (float, optional)      [%] iron loss in the MV/HV transformer - default is 0.4%
        copper_loss (float, optional)    [%] copper loss in the MV/HV transformer - default is 0.04% for trafo>50MVA, otherwise 0.2% for trafo <50MVA
   """ 

    p_nom: float|None = None
    iron_loss: float|None = None
    copper_loss: float|None = None
   
    def input_check(self):
        self.p_nom = _assign_numeric(self.p_nom, 'MV/HV transformer rated power', lo=0, standard=None, unit='MVA', custom_msg=' (no HV trafo)')
        
        # if #p_nom is not None, assign all other class attributes
        if self.p_nom is not None:
            if self.p_nom >=50:
                standard_copper_loss = 0.04
            else:
                standard_copper_loss = 0.02
        
            self.iron_loss = _assign_numeric(self.iron_loss, 'MV/HV transformer iron loss', lo=0, hi= 100, standard=0.4, unit='%')
            self.copper_loss = _assign_numeric(self.copper_loss, 'MV/HV transformer copper loss', lo=0, hi= 100, standard=standard_copper_loss, unit='%')
            # compare values to ensure copper loss are lower than iron loss.
            _compare_numeric(self.copper_loss, self.iron_loss, 'MV/HV trafo copper loss', 'MV/HV trafo iron loss', custom_msg='. Check HV trafo inputs input')       

@dataclass
class HV_Connection:
    """
    Defines HV connection parameters
    
    Attributes:
        line_len (optional, float)    [km] lenght of the HV connection cable, from HV trafo (if present) to PoC. If None provided, will assume no HV connection. Loss is ~0.1% per km of line
    """

    line_len: float|None = None

    def input_check(self):
        self.line_len = _assign_numeric(self.line_len, 'HV connection line lenght', lo=0, standard=None,  custom_msg=' (no HV connection line)')      
            
@dataclass
class Project_Curtailment:
    """
    Defines PV_project curtailment 

    Attributes:
        grid_curtail_point  (optional,str)  [-] where apply grid curtail. Available options are: (NOT case sensitive)
                                                     # None                  : No grid curtailment                                        
                                                     # 'Final MV line in'    : the grid curtailment is applied at the entry of the final MV connection line (the one that goes to HV trafo).
                                                     # 'HV trafo in'         : the grid curtailment is applied at the exit of the final MV connection line, i.e. the entry of the HV trafo (if present).
                                                     # 'HV trafo out'        : the grid curtailment is applied at HV trafo output, i.e. the entry of the HV line. Useful when PoC is right after MV/HV trafo
                                                     # 'HV line out'         : the grid curtailment is applied at HV line output. Useful when PoC is at the end of the HV line
                                                     # NOTE: 'HV trafo in' is the same physical node as the exit of the final MV line, and
                                                     #       'HV trafo out' is the same physical node as the entry of the HV line.
        grid_curtail_power (optional, float)    [MVA] grid AC power curtail for point specified
    """
    
    grid_curtail_point: str|None = None
    grid_curtail_power: float|None = None
    
    def input_check(self):
        # the value is normalised (lowercase, single spaces) and remapped through the alias table, so that user input is NOT case sensitive
        self.grid_curtail_point = _assign_str_from_list(self.grid_curtail_point, 'PV Project curtailment application point', _VALID_PROJECT_CURTAILMENT_POINT, standard = None, custom_msg=' hence no curtailment applied')
        if self.grid_curtail_point is not None:
            self.grid_curtail_power = _assign_numeric(self.grid_curtail_power, 'PV Project curtailment AC power', lo=0)

@dataclass                      
class PV_Project:
    """
    Defines parameters common to all PV fields
    
    Attributes:
        location (Location)                                  object defining project location (latitude, longitude)
        tmy (Location)                                       object defining TMY source
        pv_fields (Tuple[PV_Field])                          Tuple of PV_Field objects belonging to the project
        hv_trafo (HV_Trafo, optional)                        object defining MV/HV trafo parameters
        hv_connection (HV_Connection, optional)              object defining HV connection after MV/HV trafo
        project_curtailment (Project_Curtailment, optional)  object defining power project curtailment
        grid_availability (float, optional)                  [%] Grid availability - default is 99.7%
    """
    
    coordinates: Coordinates
    tmy: Tmy
    pv_fields: Tuple[PV_Field, ...]
    hv_trafo: HV_Trafo|None = None
    hv_connection: HV_Connection|None = None
    project_curtailment: Project_Curtailment|None = None
    grid_availability: float|None = None

    def __post_init__(self):
        """
        Replaces the optional attributes left to None with a properly initialised, empty object.

        BUGFIX: this used to be done inside input_check(), but the module level input_check(pv_project)
        function calls hv_trafo.input_check() / hv_connection.input_check() / project_curtailment.input_check()
        BEFORE pv_project.input_check(), so a project declared without HV parts crashed with an AttributeError.
        Doing it at construction time guarantees the attributes always hold a valid object.
        """

        if self.hv_trafo is None:
            self.hv_trafo = HV_Trafo()                      # no MV/HV transformer in the project

        if self.hv_connection is None:
            self.hv_connection = HV_Connection()            # no HV connection line in the project

        if self.project_curtailment is None:
            self.project_curtailment = Project_Curtailment()  # no project curtailment

    def input_check(self):
        # 0 - check that all pv_project attributes classes are good
        self.coordinates.input_check()
        self.tmy.input_check()
        for f in self.pv_fields:
            f.power.input_check(f.name)
            f.ground.input_check(f.name)
            f.mms.input_check(f.name)
            f.module.input_check(f.name)
            f.inverter.input_check(f.name)
            f.mv_trafo.input_check(f.name)
            f.field_curtailment.input_check(f.name)
            f.loss.input_check(f.name)
            f.mv_connection.input_check(f.name)
            f.input_check()
        self.hv_trafo.input_check()
        self.hv_connection.input_check()
        self.project_curtailment.input_check()
        
        # 1 - check that at least one PV_Field is provided
        if not isinstance(self.pv_fields, (tuple, list)) or len(self.pv_fields) == 0:
            raise ValueError('PV Project must contain at least one PV_Field')
        else: 
            self.pv_fields = tuple(self.pv_fields)  # in case user uses list, it is always converted to tuple
            #calculate max AC power of project as sum of max inverter power for each field
            max_project_ac_power = sum((np.ceil(f.power.p_ac_nom*1000/f.inverter.Pnom_ac)*f.inverter.P_max_ac/1000) for f in self.pv_fields)

        # 2 - check that every PV_Field carries a unique name.
        # The field name is used as a column name in the RESULTS DataFrame and as a key of the MV tree,
        # therefore two fields sharing the same name would silently overwrite each other.
        field_names = [str.lower(f.name) for f in self.pv_fields]
        duplicated_names = sorted({name for name in field_names if field_names.count(name) > 1})
        if duplicated_names:
            raise ValueError(f'PV Project: every PV field must carry a unique name. Duplicated names found: {duplicated_names}')
        
        # 3 - the optional attributes (hv_trafo, hv_connection, project_curtailment) have already been
        # initialised by __post_init__, so here they are guaranteed to hold a valid object

        # 4 - if a HV trafo is declared, check that its rating is at least 20% higher than the nominal AC power
        if self.hv_trafo.p_nom is not None:
            min_trafo_pnom = max_project_ac_power * 1.2
            if self.hv_trafo.p_nom < min_trafo_pnom:
                warnings.warn(f'PV Project: Nominal Power of HV trafo of ({self.hv_trafo.p_nom}) MVA should be at least 20% higher than the maximum PV project AC power ({max_project_ac_power} MVA), to accommodate for reactive power generation')      

        self.grid_availability = _assign_numeric(self.grid_availability, 'Grid availability', lo=0, hi= 100, standard=99.7, unit='%', custom_msg=', which is standard assumption for MV connection without step-up substation')
        

        # 5 -  checks that parameters in MV_connection are fields name or 'poc' or 'grid'
        valid_field_destination = [str.lower(f.name) for f in self.pv_fields] + _VALID_MV_TREE_EXITS
        for f in self.pv_fields:
            if str.lower(f.destination_field) not in valid_field_destination:
                raise RuntimeError (f'{f.name}: Invalid mv connection field destination. Valid destinations are {valid_field_destination}')
        
        # 6 - checks that MV connections tree is complete and that no loop arise
        # Create a temporary fields dictionary map to be used in the following two steps
        field_map = {str.lower(f.name): f for f in self.pv_fields}
        
        # --- Phase A: Tree integrity check (No loops, no dead ends) ---
        """
        this bit of codes works by:
            1) identifying a starting PV field (tree node)
            2) for that starting node, the code descends until the exit node ('grid or ' or 'poc'))
                a. if it encourters a node already visited, we have a loop!
                b. if exit node is found, without loops, that starting node is successfully connected to exit node
            3) when exit node is found, the cycle restarts from next PV_field. visited nove is resetted 
        """
        for start_field in self.pv_fields:
            visited = set()
            current_target = str.lower(start_field.destination_field)
            
            while current_target not in _VALID_MV_TREE_EXITS:
                
                # 1. Loop check
                if current_target in visited or current_target == str.lower(start_field.name):
                    raise RuntimeError(f"MV Topology Error: Loop detected! Field '{start_field.name}' "
                                       f"is stuck in an infinite cycle passing through '{current_target}'.")
                
                # 2. Dead end check (target doesn't exist in the project)
                if current_target not in field_map:
                    raise RuntimeError(f"MV Topology Error: Broken tree! Field "
                                       f"'{start_field.name}' points to a non-existent target: '{current_target}'.")
                
                # Mark as visited and continue navigating to the next parent
                visited.add(current_target)
                next_field = field_map[current_target]
                current_target = str.lower(next_field.destination_field)
                
        # --- Phase B: Physical connections construction ---
        # if loop and completeness check above works, the code assign the incoming_field attributes of each PV_field
        for f in self.pv_fields:
            dest_name = str.lower(f.destination_field)
            
            # If it has an internal destination, link this field to the destination's 'incoming_fields'
            if dest_name not in _VALID_MV_TREE_EXITS:
                target_field_obj = field_map.get(dest_name)
                if target_field_obj:
                    target_field_obj.incoming_fields.append(f)
        
        
        # 6.1 - check that exactly ONE PV field feeds the exit node of the MV tree.
        # Two fields going straight to the grid would be two separate projects, not one, and there would be
        # no single 'final MV line' where the project curtailment can be applied.
        exit_fields = [f for f in self.pv_fields if str.lower(f.destination_field) in _VALID_MV_TREE_EXITS]

        if len(exit_fields) == 0:
            raise RuntimeError('MV Topology Error: no PV field is connected to the exit of the MV tree. '
                               'Exactly one PV field must have its MV connection destination_field set to "grid" or "PoC"')
        if len(exit_fields) > 1:
            raise RuntimeError(f'MV Topology Error: {len(exit_fields)} PV fields are connected directly to the exit of the MV tree '
                               f'({[f.name for f in exit_fields]}). Only one PV field can feed the Point of Connection: '
                               'two independent exits describe two separate PV projects, which must be simulated separately')

        # store the field owning the last MV line of the tree. The engine applies the project curtailment here
        self.exit_field = exit_fields[0]

        # 6.2 - check that every PV field shares the same MV voltage.
        # All the fields are galvanically connected through the MV tree, hence they must operate at the same voltage
        mv_voltages = {f.mv_trafo.v_out for f in self.pv_fields}
        if len(mv_voltages) > 1:
            detail = {f.name: f.mv_trafo.v_out for f in self.pv_fields}
            raise ValueError(f'PV Project: all PV fields are connected to the same MV network, therefore they must share '
                             f'the same LV/MV transformer output voltage. Voltages found [kV]: {detail}')

        # 7 - check that curtailment data match the selected curtailment scenario and PV fields number
        if self.project_curtailment.grid_curtail_point is not None:
            if self.project_curtailment.grid_curtail_power > max_project_ac_power:
                warnings.warn(f'PV Project: AC power curtailment is ({self.project_curtailment.grid_curtail_power}) MVA but should be lower than the sum of maximum AC power of all inverters ({max_project_ac_power} MVA)')


class Results:
    """
    Builds and describes the RESULTS DataFrame produced by the PVyield engine.

    The RESULTS DataFrame is organised as follows:
        - index    : the loss/gain chain, in the exact order in which the losses are physically applied
        - 'Units'  : unit of measurement of each row
        - one column per PV field, holding the losses computed for that field alone
        - 'PV project' column, holding the losses of the whole project

    The chain is split into four blocks, separated by spacer rows:

        BLOCK 1 - PV FIELD
            From GHI down to the field MV cabin. Every row is computed independently for each PV field.
            Closes with 'PR fields' and 'Yield fields'.

        BLOCK 2 - PV FIELDS + MV NETWORK
            Adds PV availability and the losses of the MV connection lines that link the fields to each other
            and, ultimately, to the exit of the MV tree.
            Closes with 'PR fields + MV connection + PV avail.' and its Net Yield.

        BLOCK 3 - PROJECT CURTAILMENT AND COMMON HV PARTS
            Every project curtailment point, HV transformer, HV connection line and grid availability.
            These are common to all PV fields, therefore they are computed only for the 'PV project' column and
            left as pd.NA in the individual field columns, together with the PR and Net Yield that close the block.

        BLOCK 4 - PROJECT RESULT
            'PR fields + MV&HV connection + grid avail.' and its Net Yield.

    Sign convention:
        - every row is stored as a unitless multiplicative deviation: -0.02 means a 2% loss, +0.02 means a 2% gain
        - the Performance Ratio of a block is the running product of (1 + value) over the rows of that block
    """

    # ---------------------------------------------------------------------
    # SPACER ROWS
    # ---------------------------------------------------------------------
    # The spacers are dashed lines used to separate the blocks of the loss chain.
    #
    # NOTE: the spacers MUST be five DIFFERENT strings, otherwise pandas collapses them into a single
    # DataFrame row and the blocks lose their separation. Uniqueness is obtained by appending a growing
    # number of trailing spaces: they are invisible once printed, but they make each string unique.
    #
    # The dash count is set by _resize_spacers(), called at the end of the class definition, so that the
    # dashed line is exactly as long as the longest row name of the chain and the printed table looks tidy.
    SPACER_1 = '-'
    SPACER_2 = '-' + ' '
    SPACER_3 = '-' + '  '
    SPACER_4 = '-' + '   '
    SPACER_5 = '-' + '    '
    SPACERS = (SPACER_1, SPACER_2, SPACER_3, SPACER_4, SPACER_5)

    # ---------------------------------------------------------------------
    # ROW NAMES
    # ---------------------------------------------------------------------
    # Rows are referenced by name (never by position) so that the chain can be reordered or extended
    # without silently breaking the Performance Ratio calculation.

    GHI                 = 'GHI'
    POA_GAIN            = 'POA gain'

    # BLOCK 1 - losses from POA to the field MV cabin (the PR of a single field is the product over these rows)
    PV_ROWS = ['Horizon loss',
                  'Near shading loss',
                  'IAM loss',
                  'Soiling loss',
                  'POA gain back',
                  'Temp + Irr loss',
                  'Partial shading loss',
                  'Module quality gain',
                  'LID + LeTID loss',
                  'Mismatch front loss',
                  'Mismatch back loss',
                  'DC ohmic wiring loss',
                  'Inverter MPPT V/I loss',
                  'Inverter efficiency loss',
                  'Inverter clipping loss',
                  'Curtailment @ inverter',
                  'LV AC cable loss',
                  'LV/MV trafo loss',
                  'Auxiliaries loss',
                  'MV AC cable loss',
                  'First year degradation',
                  'Curtailment @ MV cabin',
                  ]

    PR_PV_ONLY    = 'PR PV fields'
    YIELD_PV_ONLY = 'Yield PV fields'

    # BLOCK 2 - availability and MV tree
    FIELDS_ROWS = ['PV availability',
                   'MV connection loss (field-to-field)',
                   ]

    PR_PV_MV       = 'PR fields + MV connection + PV avail.'
    YIELD_PV_MV    = 'Yield fields + MV connection + PV avail.'

    # BLOCK 3 - project curtailment and common HV parts. Computed for the project only.
    # 'Curtailment @ Final MV line in' belongs here: even though it is physically applied inside the MV tree,
    # it is a PROJECT curtailment, set through the Project_Curtailment class, and it has no meaning for a single field
    PROJECT_ROWS = ['Curtailment @ Final MV line in',
                    'Curtailment @ HV trafo in',
                    'HV trafo loss',
                    'Curtailment @ HV trafo out',
                    'HV connection loss',
                    'Curtailment @ HV line out',
                    'Grid availability',
                    ]

    PR_PROJECT      = 'PR fields + MV&HV connection + grid avail.'
    YIELD_PROJECT   = 'Yield fields + MV&HV connection + grid avail.'

    # Rows that are physically common to the whole project and therefore meaningless for a single field.
    # In the field columns these are left as pd.NA and treated as a neutral factor (0) when building the PR.
    COMMON_ROWS = list(PROJECT_ROWS)

    # Name of the project column
    PROJECT_COLUMN = 'PV project'

    @staticmethod
    def _resize_spacers():
        """
        Resizes the spacer rows so that the dashed line is exactly as long as the longest row name of the chain.

        This is run once, right after the class body is defined, so that adding or renaming a row of the loss
        chain automatically adjusts the dashes without anyone having to count characters by hand.

        The spacers keep a growing number of trailing spaces, which makes each of them a unique string:
        without it pandas would merge them into a single DataFrame row and the blocks would lose their separation.
        """

        # every row name that will appear in the index of the RESULTS DataFrame
        all_rows = ([Results.GHI, Results.POA_GAIN]
                    + Results.PV_ROWS
                    + [Results.PR_PV_ONLY, Results.YIELD_PV_ONLY]
                    + Results.FIELDS_ROWS
                    + [Results.PR_PV_MV, Results.YIELD_PV_MV]
                    + Results.PROJECT_ROWS
                    + [Results.PR_PROJECT, Results.YIELD_PROJECT])

        dashes = '-' * max(len(row) for row in all_rows)

        Results.SPACER_1 = dashes
        Results.SPACER_2 = dashes + ' '
        Results.SPACER_3 = dashes + '  '
        Results.SPACER_4 = dashes + '   '
        Results.SPACER_5 = dashes + '    '
        Results.SPACERS = (Results.SPACER_1, Results.SPACER_2, Results.SPACER_3, Results.SPACER_4, Results.SPACER_5)

    @staticmethod
    def value_columns(RESULTS):
        """
        Returns the list of columns of the RESULTS DataFrame holding numeric values.

        The 'Units' column holds strings and must be excluded from any numeric assignment, otherwise
        pandas raises a dtype error when a float is written across a whole row.

        Args:
            RESULTS (pd.DataFrame):  RESULTS DataFrame

        Returns:
            list[str]:  every column except 'Units'
        """
        return [col for col in RESULTS.columns if col != 'Units']

    @staticmethod
    def assign_results(pv_project):
        """
        Preallocates the RESULTS DataFrame and fills in every loss that is already known from the input classes.

        Args:
            pv_project (PV_Project):   PV_Project class object

        Returns:
            RESULTS (pd.DataFrame):    DataFrame indexed on the loss chain, with a 'Units' column,
                                       one column per PV field and one 'PV project' column
        """

        results_name_units = {}

        results_name_units[Results.GHI] = 'kWh/m2'
        results_name_units[Results.POA_GAIN] = '%'
        for row in Results.PV_ROWS:
            results_name_units[row] = '%'
        results_name_units[Results.SPACER_1] = '--'
        results_name_units[Results.PR_PV_ONLY] = '%'
        results_name_units[Results.YIELD_PV_ONLY] = 'kWh/kWp'

        results_name_units[Results.SPACER_2] = '--'
        for row in Results.FIELDS_ROWS:
            results_name_units[row] = '%'
        results_name_units[Results.SPACER_3] = '--'
        results_name_units[Results.PR_PV_MV] = '%'
        results_name_units[Results.YIELD_PV_MV] = 'kWh/kWp'

        results_name_units[Results.SPACER_4] = '--'
        for row in Results.PROJECT_ROWS:
            results_name_units[row] = '%'
        
        results_name_units[Results.SPACER_5] = '--'
        results_name_units[Results.PR_PROJECT] = '%'
        results_name_units[Results.YIELD_PROJECT] = 'kWh/kWp'

        RESULTS = pd.DataFrame(index=list(results_name_units.keys()))
        RESULTS['Units'] = list(results_name_units.values())

        # preallocate every result as zero, for each PV field and for the project
        column_list = [f.name for f in pv_project.pv_fields] + [Results.PROJECT_COLUMN]
        for col in column_list:
            RESULTS[col] = np.zeros(len(results_name_units))

            # spacer rows carry no value
            for spacer in Results.SPACERS:
                RESULTS.loc[spacer, col] = pd.NA

        # rows describing the common HV parts are meaningless at field level
        for f in pv_project.pv_fields:
            for row in Results.COMMON_ROWS:
                RESULTS.loc[row, f.name] = pd.NA

        # preallocate the losses that are already known from the input classes of each field
        for f in pv_project.pv_fields:
            RESULTS.loc['Soiling loss', f.name]           = -(f.loss.soiling_loss)/100
            RESULTS.loc['Module quality gain', f.name]    = -(f.loss.module_quality)/100
            RESULTS.loc['LID + LeTID loss', f.name]       = -(f.loss.lid_letid_loss)/100
            RESULTS.loc['Mismatch front loss', f.name]    = -(f.loss.mismatch_loss_front)/100
            RESULTS.loc['Mismatch back loss', f.name]     = -(f.loss.mismatch_loss_back)/100
            RESULTS.loc['Auxiliaries loss', f.name]       = -(f.loss.auxiliaries_loss)/100
            RESULTS.loc['First year degradation', f.name] = -(f.loss.first_year_degrad)/100
            RESULTS.loc['PV availability', f.name]        = f.loss.pv_availability/100-1

        # grid availability is common to the whole project
        RESULTS.loc['Grid availability', Results.PROJECT_COLUMN] = pv_project.grid_availability/100-1

        return RESULTS

    @staticmethod
    def compute_performance(RESULTS, column, glob_inc, energy, p_dc_peak):
        """
        Computes the three Performance Ratios and the three Net Yields of one column of the RESULTS DataFrame.

        The Performance Ratio is computed the same way PVsyst does, as the ratio between the specific energy
        actually delivered and the irradiation received by the plane of array:

            Net Yield [kWh/kWp] = energy / peak DC power
            PR                  = Net Yield / GlobInc

        where GlobInc is the global irradiation on the collector plane, front side only, before the horizon
        and the near shading losses. That is exactly the reference PVsyst uses, which makes the two Performance
        Ratios comparable one to one.

        NOTE: the PR is NOT the product of the loss rows. The rows describe where the energy goes, but two
        coefficients convert the optical chain into the electrical one, the STC efficiency and the bifaciality
        factor, and neither of them is a loss. PVsyst applies both at the 'PV conversion' step of its loss
        diagram, and computes its PR separately as Yf/Yr, which is what is reproduced here. As a consequence
        the bifaciality raises the PR, because it enters the delivered energy but not the reference irradiation.

        Args:
            RESULTS (pd.DataFrame):  RESULTS DataFrame
            column (str):            name of the column to be processed (a PV field name, or 'PV project')
            glob_inc (float):        [kWh/m2] global irradiation on the collector plane, front side, single shed
            energy (dict):           [Wh] energy delivered at the end of each block, with keys 1, 2 and 3.
                                     Block 3 only exists for the project column and can be left out
            p_dc_peak (float):       [Wp] peak DC power the Net Yield is normalised on

        Returns:
            RESULTS (pd.DataFrame):  the same DataFrame, with the PR and Net Yield rows filled in
        """

        def block_result(block, pr_row, yield_row):
            """Net Yield in kWh/kWp and PR in %, from the energy delivered at the end of one block."""
            if block not in energy or not p_dc_peak or not glob_inc:
                RESULTS.loc[pr_row, column] = pd.NA
                RESULTS.loc[yield_row, column] = pd.NA
                return

            net_yield = energy[block]/p_dc_peak            # [Wh/Wp] = [kWh/kWp]
            RESULTS.loc[yield_row, column] = net_yield
            RESULTS.loc[pr_row, column] = net_yield/glob_inc*100

        # BLOCK 1 - single PV field, up to the field MV cabin
        block_result(1, Results.PR_PV_ONLY, Results.YIELD_PV_ONLY)

        # BLOCK 2 - adds availability and the MV connection tree
        block_result(2, Results.PR_PV_MV, Results.YIELD_PV_MV)

        # BLOCK 3 - adds the project curtailment and the common HV parts. These belong to the shared part of
        # the plant, so they only make sense for the project: in the field columns they are left as pd.NA
        block_result(3, Results.PR_PROJECT, Results.YIELD_PROJECT)

        return RESULTS

    @staticmethod
    def print_results(RESULTS):
        """
        Prints the loss chain of the PROJECT on the console.

        Only three columns are shown: the name of the loss, its unit of measurement and the value computed
        for the whole project. The columns of the single PV fields are left out on purpose: they would make
        the table too wide to be read on a console, and they are available in full in the Excel report.

        Args:
            RESULTS (pd.DataFrame):  RESULTS DataFrame, already formatted through format_results()
        """

        # keep the loss name (the index), the unit and the project column only
        results_to_print = RESULTS[['Units', Results.PROJECT_COLUMN]]

        # DISPLAY OPTIONS
        # The three dots that used to appear between the units and the values were pandas TRUNCATING THE COLUMNS:
        # when a DataFrame holds more columns than 'display.max_columns', pandas replaces the middle ones with '...'.
        # Spyder and IPython set that option to a small number of their own, so the more PV fields in the project,
        # the sooner the table got truncated. The options below disable every truncation, on rows and columns alike.
        # option_context restores the previous settings on exit, so the user session is left untouched.
        with pd.option_context('display.float_format', '{:+,.2f}'.format,
                               'display.max_rows', None,            # never truncate the rows of the loss chain
                               'display.max_columns', None,         # never truncate the columns (this removes the '...')
                               'display.width', 1000,               # do not wrap on the console width
                               'display.expand_frame_repr', False,  # keep every column on the same line, no block wrapping
                               ):
            print(results_to_print)

    @staticmethod
    def format_results(RESULTS):
        """
        Converts the loss rows from unitless to % and rounds the whole DataFrame for printing.

        Args:
            RESULTS (pd.DataFrame):  RESULTS DataFrame holding unitless values

        Returns:
            RESULTS (pd.DataFrame):  RESULTS DataFrame with loss rows expressed in %
        """

        # every row expressed in % is scaled by 100. GHI, Net Yields and spacers are left untouched
        percent_rows = [Results.POA_GAIN] + Results.PV_ROWS + Results.FIELDS_ROWS + Results.PROJECT_ROWS

        for col in Results.value_columns(RESULTS):
            RESULTS.loc[percent_rows, col] = pd.to_numeric(RESULTS.loc[percent_rows, col], errors='coerce') * 100
            RESULTS[col] = pd.to_numeric(RESULTS[col], errors='coerce').round(2)

        return RESULTS




# Resize the spacer rows on the longest row name of the loss chain.
# It is called here, right after the class definition, because a staticmethod cannot be invoked
# from inside the class body while the class itself is still being built.
Results._resize_spacers()


# ==========================================================================
# SUPPORT FUNCTIONS
# ==========================================================================


class WarningLog:
    """
    Collects every warning raised while the simulation runs, so that they can be reported in the Excel file.

    The Python warnings module prints its messages to the console and forgets them. This class installs a hook
    that keeps printing them exactly as before, but also stores each message in a list, which excel_generator()
    then writes into a dedicated sheet of the report.

    Usage:
        WarningLog.start()                     # called once, at the very beginning of the launcher
        ... build the classes, run the simulation, write the report ...
        messages = WarningLog.get_messages()   # list of the warnings raised, in the order they appeared
        WarningLog.stop()                      # called once, at the very end of the launcher

    start() and stop() must always come in pairs. While the log is running, the warnings module is left in a
    modified state (a custom showwarning hook and the 'always' filter), which would otherwise survive the run
    and affect any other code executed later in the same console session, as happens in Spyder or Jupyter
    where the interpreter stays alive between runs.

    Attributes:
        messages (list[dict]):  one entry per warning, holding its category, its message and the name of the
                                source file and line that raised it
    """

    messages = []
    _active = False
    _previous_showwarning = None
    _previous_filters = None

    @staticmethod
    def start():
        """
        Starts collecting the warnings, and clears anything collected by a previous run.

        The default warning filter of Python only shows a given warning ONCE per source line, so the second
        PV field raising the very same message at the very same line would be silently dropped. The filter is
        therefore set to 'always', which is also what the user expects to see when several PV fields share the
        same missing input.
        """

        WarningLog.messages = []

        if WarningLog._active:
            return   # the hook is already installed, nothing else to do

        # remember how the warnings module was configured, so that stop() can put it back as it was
        WarningLog._previous_showwarning = warnings.showwarning
        WarningLog._previous_filters = warnings.filters[:]

        # show every warning, including the repeated ones coming from different PV fields.
        # Without this, the default Python filter would show a given message only ONCE per source line,
        # and the same warning raised by a second PV field would be silently dropped
        warnings.simplefilter('always')

        def _collect(message, category, filename, lineno, file=None, line=None):
            WarningLog.messages.append({'category': category.__name__,
                                        'message': str(message),
                                        'source': f'{Path(filename).name}, line {lineno}',
                                        })
            # keep the original behaviour, so the warnings still reach the console
            WarningLog._previous_showwarning(message, category, filename, lineno, file, line)

        warnings.showwarning = _collect
        WarningLog._active = True

    @staticmethod
    def stop():
        """
        Stops collecting the warnings and restores the standard behaviour of the warnings module.

        Both the showwarning hook and the warning filters are put back as they were before start() was called,
        so that the rest of the console session is left exactly as the log had found it.

        The messages already collected are kept, so they can still be written to the Excel report.
        """

        if not WarningLog._active:
            return   # nothing to restore

        warnings.showwarning = WarningLog._previous_showwarning

        # restore the filters, otherwise the simplefilter('always') set by start() would survive the run
        if WarningLog._previous_filters is not None:
            warnings.filters[:] = WarningLog._previous_filters
            warnings._filters_mutated()   # tells the warnings module to drop its internal cache

        WarningLog._active = False

    @staticmethod
    def get_messages():
        """
        Returns the warnings collected so far.

        Returns:
            list[dict]:  one entry per warning, each holding the keys 'category', 'message' and 'source'
        """

        return list(WarningLog.messages)


def custom_formatwarning(message, category, filename, lineno, line=None):
    """
    Custom formatter for the warnings module to output a cleaner message.
    
    Args:
        message (Warning or str): The warning message to be displayed.
        category (type): The warning category/class (e.g., UserWarning, DeprecationWarning).
        filename (str): The name of the file where the warning was triggered.
        lineno (int): The line number in the source code where the warning occurred.
        line (str, optional): The line of source code that caused the warning. Defaults to None.
        
    Returns:
        str: A formatted string containing only the warning category and the message,
             without the file path, line number, or repeated source code.
    """
    # Format the string to return only "Category: Message"
    return f"{category.__name__}: {message}\n\n"


def _assign_numeric(value: any,
                   name: str,
                   lo: float|None = None,
                   hi: float|None = None,
                   standard: float|int|None  = _EMPTY,
                   unit: str = '',
                   int_only: bool = False,
                   custom_msg: str = '',
                   ):
    """
    Validates that given numeric value is int or float and that it falls within a specified range. If valid, return the value. If not valid, returns the fallsback value (if present)

    Args:
        value (any):                The numeric value to be checked.
        name (str):                 The name of the variable (used for formatting error/warning messages).
        lo (float, optional):       The lower boundary. If provided, `value` must be >= `lo`. Default is None.
        hi (float, optional):       The upper boundary. If provided, `value` must be <= `hi`. Default is None.
        standard (float, optional): The fallback value to return if the check fails.  If provided, a warning is issued. If not provided, the fallback value is returned.
        unit (str, optional):       The unit of measurement for the value (e.g., 'MW' or 'deg'). Only required if standard is defined
        int_only (bool, optional):  If True, function will accept only int type. If false, function will accept float and int types
        custom_msg (str, optional): Custom message to be appended in warning message, to give details about standard value used
        
    Returns:
        float or int: The original `value` if it passes the checks, or the `standard` value if it fails.
    """
    
    if int_only:
        valid_types = int
        type_str = 'integer'
    else:
        valid_types = (int, float)
        type_str = 'int or float'
    
    # peculiar case: value is None. If Value is None, the function will crash since None cannot be used with operators (<>), hence, None case is treated aside
    if value is None:
        if standard is None:
            return value
        elif standard is _EMPTY: # standard is not assigned in when calling the function
            raise TypeError(f'{name} must be a number {type_str}. Please check inputs and provide a correct value')
        else: # standard is a legit value
            warnings.warn(f'Invalid or missing {name}. Value is None. The simulation assumed a value of {standard} {unit}{custom_msg}')
            return standard
    
    #TYPE CHECK - check that value is correct type. Check on bool is needed because bool values are subclass of int (bool is 0 or 1 indeed)
    if not isinstance(value, valid_types) or isinstance(value, bool):
        if standard is not _EMPTY:
            warnings.warn(f'WARNING: {name} must be a numeric {type_str}, not {type(value).__name__}. The simulation assumed a value of {standard} {unit}{custom_msg}')
            return standard
        else:
            raise TypeError(f'{name} must be a number {type_str}, not {type(value).__name__}. Please check inputs and provide a correct value')


    # CASE A: both boundaries provided. Function check that value is between hi and lo
    if lo is not None and hi is not None:
        if not (lo <= value <= hi):
            if standard is not _EMPTY:
                warnings.warn(f'Invalid or missing {name}. Value must be between {lo} {unit} and {hi} {unit}. The simulation assumed a value of {standard} {unit}{custom_msg}')
                return standard
            else:
                raise ValueError(f'Invalid or missing {name}. Value must be between {lo} {unit} and {hi} {unit}. Please check inputs and provide a correct value')
                 
    # CASE B: only lo is provided. Function check that value is higher than lo
    elif lo is not None:
        if value < lo:
            if standard is not _EMPTY:
                warnings.warn(f'Invalid or missing {name}. Value must be higher than {lo} {unit}. The simulation assumed a value of {standard} {unit}{custom_msg}')
                return standard
            else:
                raise ValueError(f'Invalid or missing {name}. Value must be higher than {lo} {unit}. Please check inputs and provide a correct value')
            
            
    # CASE C: only hi is provided. Function check that value is lower than hi
    elif hi is not None:
        if value > hi:
            if standard is not _EMPTY:
                warnings.warn(f'Invalid or missing {name}. Value must be lower than {hi} {unit}. The simulation assumed a value of {standard} {unit}{custom_msg}')
                return standard
            else:
                raise ValueError(f'Invalid or missing {name}. Value must be lower than {hi} {unit}. Please check inputs and provide a correct value')
                
    # CASE 2.D : neither lo or hi are provided, or if value is valid and within boundaries
    elif hi is None and lo is None:
            return value
    
    return value


def _assign_str (value: any,
                 name: str,
                 standard: str|None = _EMPTY,
                 custom_msg: str = '',
                 ):
    """
    Validates that given value is of type string. If valid, return the value. If not valid, returns the fallsback value (if present)

    Args:
        value (str): The value to be checked.
        name (str):  The name of the variable (used for formatting error/warning messages).
        standard (str, optional): The fallback value to return if the type check fails. If provided, a warning is issued. If None, an exception is raised. Default is None.
        custom_msg (str, optional): Custom message to be appended in warning message, to give details about standard value used

    Returns:
        str: The original `value` if it is a string, or the `standard` value if the check fails.
    """
    
    if not isinstance(value, str):
        if standard is not _EMPTY:
            warnings.warn(f'Invalid or missing {name}. This must be a string, but you provided a {type(value).__name__}. The simulation assumed {standard}{custom_msg}')
            return standard
        else:
            raise ValueError(f'Invalid or missing {name}. This must be a string. Please check inputs and provide a correct value')
    return value


def _assign_numeric_from_list(value: any,
                               name: str,
                               valid_list: list[int|float],
                               standard: int|float|None = _EMPTY,
                               unit: str = '',
                               custom_msg: str = '',
                               ):
    """
    Validates that a float value belongs to a specified list of allowed options. If valid, return the value. If not valid, returns the fallsback value (if present)

    Args:
        value (float):               The float value to be checked.
        name (str):                  The name of the variable (used for formatting error/warning messages).
        valid_list (List[float]):    A list containing all the acceptable float values.
        standard (str, optional):    The fallback value to return if the check fails.  If provided, a warning is issued. If None, an exception is raised. Default is None.
        unit (str, optional):        The unit of measurement for the value (e.g., 'MW' or 'deg'). Only required if standard is defined
        custom_msg (str, optional):  Custom message to be appended in warning message, to give details about standard value used

    Returns:
        float: The original `value` if it is valid, or the `standard` value if it fails the check.
    """
    
    if not isinstance (value, (int, float)) or value not in valid_list:
        if standard is not _EMPTY:
            warnings.warn(f'Invalid or missing {name}. Valid values are: {valid_list}. The simulation assumed a value of {standard} {unit}{custom_msg}')
            return standard
        else:
            raise ValueError(f'Invalid or missing {name}. Valid values are: {valid_list}. Please check inputs and provide a correct value')
    return value  


def _assign_numeric_array_like(value: any,
                               name: str,
                               exact_size: int|Tuple[int, ...] | None = None,
                               lo: float|int|None = None,
                               hi: float| int|None = None,
                               standard: float|int = _EMPTY,
                               unit: str = '',
                               int_only: bool = False,
                               custom_msg: str = '',
                               ):
    """
    Validates that an input is a numeric array-like object (list, tuple, np.ndarray). 
    It can also check for shape constraints, numerical boundaries, and integer types.
    If valid, returns the original value. If not valid, returns the fallback value.

    Args:
        value (Any):                                The array-like value to be checked.
        name (str):                                 The name of the variable (used for formatting error messages).
        exact_size (int | tuple | None, optional):  Checks the n-dimensional shape (e.g., 3 for 1D, (2, 2) for 2D).
        lo (float | int | None, optional):          Lower bound. If provided, all elements must be >= lo.
        hi (float | int | None, optional):          Upper bound. If provided, all elements must be <= hi.
        standard (Any, optional):                   The fallback value to return if the check fails. Default is _EMPTY.
        unit (str, optional):                       Unit of measurement for the warning message.
        int_only (bool, optional):                  If True, ensures that all elements are integers. Default is False.
        custom_msg (str, optional):                 Custom message to be appended to the warning message.

    Returns:
        The original value if it is valid, or the standard value if it fails the check.
    """
    
    # 1. Prepare dynamic error messages based on provided constraints
    size_err = f' of shape {exact_size}' if exact_size is not None else ''
    type_err = 'integer ' if int_only else 'numeric '
    
    bounds_err = ''
    if lo is not None and hi is not None:
        bounds_err = f' with values between {lo} and {hi}'
    elif lo is not None:
        bounds_err = f' with values >= {lo}'
    elif hi is not None:
        bounds_err = f' with values <= {hi}'
    
    # 2. Verify that it is an accepted array-like object
    if isinstance(value, (tuple, list, np.ndarray)):
        try:
            # Temporarily convert to a numpy array for vectorised checks
            arr = np.array(value)
            
            # Check that it is not empty
            if arr.size > 0:
                
                # Check data type (integer only or any valid number)
                is_valid_type = np.issubdtype(arr.dtype, np.integer) if int_only else np.issubdtype(arr.dtype, np.number)
                
                if is_valid_type:
                    # Check numerical boundaries (if any value violates the bounds, it fails)
                    bounds_ok = True
                    if lo is not None and np.any(arr < lo):
                        bounds_ok = False
                    if hi is not None and np.any(arr > hi):
                        bounds_ok = False
                        
                    if bounds_ok:
                        # Check the n-dimensional shape if requested
                        if exact_size is not None:
                            target_shape = (exact_size,) if isinstance(exact_size, int) else tuple(exact_size)
                            if arr.shape == target_shape:
                                return value
                        else:
                            return value  # All checks passed, return original value
                            
        except Exception:
            # If conversion or checks fail, fall through to the fallback logic
            pass
                
    # 3. Validation failed: return standard if provided, else raise Exception
    if standard is not _EMPTY:
        warnings.warn(f'Invalid or missing {name}. The simulation assumed a value of {standard} {unit}{custom_msg}')
        return standard
    else:
        raise ValueError(
            f'Invalid or missing {name}. The value should be a {type_err}array-like '
            f'(tuple, list, or ndarray){size_err}{bounds_err}. Please check inputs and provide a correct value.'
        )


def _assign_str_from_list(value: any,
                           name: str,
                           valid_list: List[str|None],
                           standard: str|None = _EMPTY,
                           custom_msg: str = '',
                           ):
    """
    Validates that a string value belongs to a specified list of allowed options. If valid, return the value. If not valid, returns the fallsback value (if present)

    Args:
        value (any):                The string value to be checked.
        name (str):                 The name of the variable (used for formatting error/warning messages).
        valid_list (Tuple(str)):    A list containing all the acceptable string values.
        standard (str, optional):   The fallback value to return if the check fails.  If provided, a warning is issued. If None, an exception is raised. Default is None.
        custom_msg (str, optional): Custom message to be appended in warning message, to give details about standard value used

    Returns:
        str: The original `value` if it is valid, or the `standard` value if it fails the check.
    """
    
    # if value is None and in valid_list, None is returned. this is necessary because None will not be recognised as a string
    if value is None and value in valid_list:
        return value
    
    elif not isinstance(value, str):
        if standard is not _EMPTY:
            warnings.warn(f'Invalid {name} type. Valid types are: {valid_list}. The simulation assumed: {standard}{custom_msg}')
            return standard
        else:
            raise ValueError(f'Invalid {name} type. Valid values are: {valid_list}. Please check inputs and provide a correct value')
    
    elif str.lower(value) not in valid_list:
        if standard is not _EMPTY:
            warnings.warn(f'Invalid or missing {name} type. Valid types are: {valid_list}. The simulation assumed: {standard}{custom_msg}')
            return standard
        else:
            raise ValueError(f'Invalid or missing {name} type. Valid values are: {valid_list}. Please check inputs and provide a correct value')
    
    return str.lower(value) 

    

def _assign_bool (value: any,
                  name: str,
                  standard: bool|None = _EMPTY,
                  custom_msg: str = '',
                  ):
    """
    Validates that given value is of type boolean. If valid, return the value. If not valid, returns the fallsback value (if present)

    Args:
        value (Any):                The value to be checked.
        name (str):                 The name of the variable (used for formatting error/warning messages).
        standard (bool, optional):  The fallback value to return if the type check fails. If provided, a warning is issued. If None, an exception is raised. Default is None.
        custom_msg (str, optional): Custom message to be appended in warning message, to give details about standard value used

    Returns:
        bool: The original `value` if it is a boolean, or the `standard` value if the check fails.
    """
    
    
    if not isinstance(value, bool):
        if standard is not _EMPTY:
            warnings.warn(f'Invalid {name}. This must be True or False. The simulation assumed a value of {standard}{custom_msg}')
            return standard
        else:
            raise ValueError(f'Invalid {name}. This must be True or False. Please check inputs and provide a correct value')
    return value    


def _compare_numeric (lo: float,
                     hi: float,
                     lo_name: str,
                     hi_name: str,
                     custom_msg: str = '',
                     ):

    """
    Compares two numbers, lo and hi, and checks that lo<hi

    Args:
        lo (float):       lower value
        hi (float):       higher value
        lo_name (str):    name of lower value
        hi_name (str):    name of higher value
        custom_msg (str, optional): custom message to be appended in warning message

    Returns:
        ValueError in case lo>hi
    """
    
    if lo is None or hi is None:
        raise ValueError(f'Cannot compare {lo_name} and {hi_name}: one or both values are None')
    
    if lo>hi:
        raise ValueError(f'{lo_name} must be lower than {hi_name}{custom_msg}')



def _pan_import(file_path,
               PAN_parameters: list[str]|None = None,
               pvlib_remap: bool = True,
               ):
    """
    Reads a .PAN text file and extracts variables that compare in the file as `variable_name = value`.
    Variable extracted are the ones included in PAN_parameter list (see accepted parameters list below)
    Infers the type of each variable (integer, float, boolean, string, list).
    Converts the value accordingly and returns either individual variables or a dictionary.
     
    
    Args:
        file_path : str or path object     Name or path of a .pan file
        PAN_parameters : list, optional    List of variable names to extract from the .PAN file. list must contains names that match the ones in PAN file. If None (default), a standard comprehensive list is used.
        pvlib_remap: bool, optional        If True, remap dictionary key to match PVlib library convention. If False, PAN file names will be used (currently, no naming convention exist for PVlib. remap simply remap the variable to work with PVyield convention)
    
    
    Accepted PAN_parameters are:
    -------------------
    Manufacturer (str)        [-] PV module manufacturer
    Model (str)               [-] PV module model
    PNom (float)              [Wp] Module DC power at STC
    Height (float)            [m] PV module height
    Width (float)             [m] PV module width
    Depth (float)             [m] PV module thickness
    Weight (float)            [kg] PV module weight
    NCelS (int)               [-] Number of cells connected in series. Assume full-cut cells equivalent (for half-cut, divide actual number by 2)
    NCelP (int)               [-] Number of cells in parallel. Typically 2 for half-cut and 1 for full-cut
    NDiode (int)              [-] Number of bypass diodes, typically 3
    Isc (float)               [A] Light-generated current (photocurrent) at STC / short circuit current
    Voc (float)               [V] Open circuit voltage at STC
    Imp (float)               [A] MPP current at STC
    Vmp (float)               [V] MPP voltage at STC
    muISC (float)             [mAmps/Celsius] Short-circuit current temp coefficient. On PAN is 'muIsc' (e.g. 6.57 for Trina NEG 685)
    Gamma (float)             [-] Diode ideality factor. 1.02 is standard PVsyst model
    muGamma (float)           [1/K] Temperature coefficient for the diode ideality factor
    RShunt (int)              [ohm] Shunt resistance at STC. If None, inferred from other inputs
    Rp_0 (int)                [ohm] Shunt resistance at zero irradiance. If None, inferred
    Rp_Exp (float)            [-] Single diode exponential parameter. Default a 5.5
    RSerie (float)            [ohm] Series resistance at STC. If None, inferred
    Absorb (float)            [-] Module absorptivity (1 - reflectivity)
    BifacialityFactor (float) [-] Bifaciality factor
    iam_theta (list[int])     [degrees] Array of angles at which the IAM is known
    iam_loss (list[float])    [-] IAM values for each angle in iam_theta
    VRevDiode (float)         [V] Bypass diode voltage drop. Default a 0.7 (for silica PV modules)
    RDiode (float)            [ohm] Bypass diode resistance
    SubModuleLayout (str)     [str] Either 'TwinHalfCells' or 'FullCells'
    PNomTolLow (float)        [Wp] power bin lower threeshold
    PNomTolUp (float)         [Wp] power bin upper threeshold
    FrontSurface (str)        [-] coating of front glass


    Returns:
        MODULE : dict       Dictionary containing PAN file variables in the order specified by keep_variable, with their values and datatype. The keys of the dictionary are the parameters names
    """
    
    # Mapping from PAN variable names to PVlib naming conventions
    pvlib_map = {'Manufacturer'      : 'manufacturer',
                 'Model'             : 'model',
                 'PNom'              : 'Pnom',               # in PVlib is Pdc_0, however I prefer Pnom
                 'Height'            : 'height',
                 'Width'             : 'width',
                 'Depth'             : 'thickness',          # not used or defined in PVlib yet
                 'Weight'            :  'weight',            # not used or defined in PVlib yet
                 'NCelS'             : 'cells_in_series',
                 'NCelP'             : 'cells_in_parallel',  # not used or defined in PVlib yet
                 'NDiode'            : 'n_diodes',
                 'Isc'               : 'I_L_ref',
                 'Voc'               : 'Voc',
                 'Imp'               : 'imp_ref',
                 'Vmp'               : 'vmp_ref',
                 'muISC'             : 'alpha_sc',
                 'Gamma'             : 'gamma_ref',
                 'muGamma'           : 'mu_gamma',
                 'RShunt'            : 'R_sh_ref',
                 'Rp_0'              : 'R_sh_0',
                 'Rp_Exp'            : 'R_sh_exp',
                 'RSerie'            : 'R_s',
                 'Absorb'            : 'alpha_absorption',
                 'BifacialityFactor' : 'bifaciality',
                 'iam_theta'         : 'theta_ref',
                 'iam_loss'          : 'iam_ref',
                 'VRevDiode'         : 'V_bypass_diode',  # not used or defined in PVlib yet
                 'RDiode'            : 'R_bypass_diode',  # not used or defined in PVlib yet     
                 'SubModuleLayout'   : 'substring_type',  # not used or defined in PVlib yet
                 'PNomTolLow'        : 'Pnom_low',        # not used or defined in PVlib yet
                 'PNomTolUp'         : 'Pnom_high',       # not used or defined in PVlib yet
                 'FrontSurface'      : 'glass_type'       # not used or defined in PVlib yet
                 }
    
    # default variable list if none is provided
    if PAN_parameters is None:
        PAN_parameters = ['Manufacturer', 'Model', 'PNom', 'Height', 'Width', 'Weight', 'Depth', 'NCelS', 'NCelP', 'NDiode', 'Isc', 'Voc', 'Imp', 'Vmp', 'muISC', 'Gamma', 'muGamma', 'RShunt', 'Rp_0', 'Rp_Exp', 'RSerie', 'Absorb', 'BifacialityFactor', 'iam_theta', 'iam_loss', 'SubModuleLayout', 'VRevDiode', 'RDiode', 'PNomTolLow', 'PNomTolUp', 'FrontSurface']

    
    
    with open(file_path, 'r') as file:
        MODULE = {}  # creates Dictionary where PAN will be extracted
        for line in file:
            # check if the ''End of PVObject pvIAMIAMProfile'' is in the line.
            # If so, stores the iam angle and iam loss in separated variables
            if 'End of TCubicProfile' in line:
                MODULE['iam_theta'] = np.zeros(shape = MODULE.get('NPtsEff',0))
                MODULE['iam_loss'] = np.zeros(shape = MODULE.get('NPtsEff',0))
                for n in np.arange(1, MODULE['NPtsEff']+1):
                    key = 'Point_' + str(n)
                    MODULE['iam_theta'][n-1] = MODULE.get(key)[0]
                    MODULE['iam_loss'][n-1] = MODULE.get(key)[1]
            
            elif '=' in line:  # Look for the presence of '='
                key, value = line.split('=', 1)  # Split into key and value
                value_cleaned = value.strip()  # Clean whitespace
                var_value = _infer_and_convert(value_cleaned)  # Get converted value
                MODULE[key.strip()] = var_value  # Store only the converted valu

    # remove all unused variables from MODULE dict
    MODULE = {key: value for key, value in MODULE.items() if key in PAN_parameters}
    
    # Remap keys to PVlib naming conventions if requested
    if pvlib_remap:
        MODULE = {pvlib_map.get(k, k): v for k, v in MODULE.items()}

    
    return MODULE


def _ond_import(file_path,
               OND_parameters: list[str]|None = None,
               pvlib_remap: bool = True,
               ):
    """
    Reads a .OND text file and extracts variables that compare in the file as `variable_name = value`.
    Variable extracted are the ones included in OND_parameter list (see accepted parameters list below)
    Infers the type of each variable (integer, float, boolean, string, list).
    Converts the value accordingly and returns either individual variables or a dictionary.
     
    
    Args:
        file_path : str or path object     Name or path of a .ond file
        OND_parameters : list, optional    List of variable names to extract from the .OND file. list must contains names that match the ones in OND file. If None (default), a standard comprehensive list is used.
        pvlib_remap: bool, optional        If True, remap dictionary key to match PVlib library convention. If False, OND file names will be used
    
    
    Accepted OND_parameters are:
    -------------------    
    Manufacturer (str)              [-] PV module manufacturer
    Model (str)                     [-] PV module model
    PNomConv (float)                [kW] nominal AC output power for inverter
    PMaxOUT (float)                 [kW] max AC output power for inverter
    VMppMin (float)                 [V] min MPPT voltage to allow MPPT to work 
    VMPPMax (float)                 [V] max MPPT voltage to allow MPPT to work
    IMaxDC (float)                  [A] max DC input current per whole inverter
    NbInputs (int)                  [-] number of maximum string per inverter
    NbMPPT (int)                    [-] number of inverter MPPT
    VOutConv (float)                [V] inverter output voltage
    TPNom (float)                   [°C] outdoor temperature at which the inverter produces the nominal power - default: 40
    TPMax (float)                   [°C] outdoor temperature at which the inverter produces the max power - default: 30
    VNomEff1 (float)                [V] inverter output voltage at which the efficiency is known for p_in_VNomEff1 and p_out_VNomEff1
    VNomEff2 (float)                [V] inverter output voltage at which the efficiency is known for p_in_VNomEff2 and p_out_VNomEff2
    VNomEff3 (float)                [V] inverter output voltage at which the efficiency is known for p_in_VNomEff3 and p_out_VNomEff3
    EfficEuroV1 (float)             [%] european efficiency for Voltage VNomEff1
    EfficEuroV2 (float)             [%] european efficiency for Voltage VNomEff2
    EfficEuroV3 (float)             [%] european efficiency for Voltage VNomEff3
    p_in_VNomEff1 (list[float])     [kWp] inverter input DC power producing an putput power of p_out_VNomEff1 at VNomEff1. Expressed as a list, each item of the list correspond to respective item in p_out_VNomEff1
    p_in_VNomEff2 (list[float])     [kWp] inverter input DC power producing an putput power of p_out_VNomEff2 at VNomEff2. Expressed as a list, each item of the list correspond to respective item in p_out_VNomEff2
    p_in_VNomEff3 (list[float])     [kWp] inverter input DC power producing an putput power of p_out_VNomEff3 at VNomEff3. Expressed as a list, each item of the list correspond to respective item in p_out_VNomEff3
    p_out_VNomEff1 (list[float])    [kWac] inverter output AC power resulting from input power of  p_in_VNomEff1 at VNomEff1. Expressed as a list, each item of the list correspond to respective item in p_in_VNomEff1
    p_out_VNomEff2 (list[float])    [kWac] inverter output AC power resulting from input power of  p_in_VNomEff2 at VNomEff2. Expressed as a list, each item of the list correspond to respective item in p_in_VNomEff2
    p_out_VNomEff3 (list[float]     [kWac] inverter output AC power resulting from input power of  p_in_VNomEff3 at VNomEff3. Expressed as a list, each item of the list correspond to respective item in p_in_VNomEff3
    TPLim1 (float)                  [°C] outdoor temperature that derates inverter AC power to PLim1
    TPLimAbs (float)                [°C] outdoor temperature that derates inverter AC power to PLimAbs
    PLim1 (float)                   [kW] inverter derated max AC power output when outdoor temperature is equal to TPLim1
    PLimAbs (float)                 [kW] inverter derated max AC power output when outdoor temperature is equal to TPLimAbs
    TanPhiMin (float)               [-] minimum power factor cos(phi)         
    TanPhiMax (float)               [-] maximum power factor cos(phi) 
    Night_Loss (float)              [W] night consumption of inverter


    Returns
        INVERTER : dict        Dictionary containing OND file variables in the order specified by keep_variable, with their values and datatype.The keys of the dictionary are the parameters names
    """

    # Mapping from PAN variable names to PVylib naming conventions
    # currently, there are no inverter namin convention for inverter variables in PVlib. 
    # This section just remap variable to match the ones used in PVyield
    pvlib_map = {'Manufacturer'      : 'manufacturer',
                 'Model'             : 'model',
                 'PNomConv'          : 'Pnom_ac',
                 'PMaxOUT'           : 'P_max_ac',
                 'VMppMin'           : 'v_min_mppt',
                 'VMPPMax'           : 'v_max_mppt',
                 'IMaxDC'            : 'i_dc_max',
                 'NbInputs'          : 'n_string_inputs',
                 'NbMPPT'            : 'n_mppt',
                 'VOutConv'          : 'v_out',
                 'TPNom'             : 'T_p_nom',
                 'TPMax'             : 'T_p_max',
                 'VNomEff1'          : 'v_eff1',
                 'VNomEff2'          : 'v_eff2',
                 'VNomEff3'          : 'v_eff3',
                 'EfficEuroV1'       : 'euro_v_eff1',
                 'EfficEuroV2'       : 'euro_v_eff2',
                 'EfficEuroV3'       : 'euro_v_eff3',
                 'p_in_VNomEff1'     : 'p_in_v_eff1',
                 'p_in_VNomEff2'     : 'p_in_v_eff2',
                 'p_in_VNomEff3'     : 'p_in_v_eff3',
                 'p_out_VNomEff1'      : 'p_out_v_eff1',
                 'p_out_VNomEff2'      : 'p_out_v_eff2',
                 'p_out_VNomEff3'      : 'p_out_v_eff3',
                 'TPLim1'            : 'T_lim_1',
                 'TPLimAbs'          : 'T_lim_2',
                 'PLim1'             : 'p_out_T_lim_1',
                 'PLimAbs'           : 'p_out_T_lim_2',
                 'TanPhiMin'         : 'cosphi_min',
                 'TanPhiMax'         : 'cosphi_max',
                 'Night_Loss'        : 'P_night',
                 } 
        
    if OND_parameters is None:
        OND_parameters = ['Manufacturer', 'Model', 'PNomConv', 'PMaxOUT', 'VMppMin', 'VMPPMax', 'IMaxDC', 'NbInputs', 'NbMPPT','VOutConv','TPNom', 'TPMax', 'VNomEff1', 'VNomEff2', 'VNomEff3', 'EfficEuroV1', 'EfficEuroV2', 'EfficEuroV3', 'p_in_VNomEff1', 'p_in_VNomEff2', 'p_in_VNomEff3', 'p_out_VNomEff1', 'p_out_VNomEff2', 'p_out_VNomEff3','TPLim1', 'TPLimAbs', 'PLim1', 'PLimAbs', 'TanPhiMin', 'TanPhiMax', 'Night_Loss']
        

    # first parse of the file
    # PROBLEM: Inverter efficiency variable have all same file name, so the parser will overwrite the first variable with the latter in the file having same name
    # SOLUTION: For each line, check if a section where inverter efficiency is start of finished
    with open(file_path, 'r') as file:
        INVERTER = {}
        V_eff_section = 0
        
        for line in file:
            # check if the VNomEff is in the line.
            # If so, stores the voltage at which the efficiency is computed in separated variables
            if 'VNomEff' in line:
                key, VNomEff_list = line.split('=', 1)  # Split into key and value
                VNomEff_list = VNomEff_list.strip().split(',')
                for i, value in enumerate(VNomEff_list, start=1):
                    variable_name = 'VNomEff'+ str(i)
                    if value.strip():      # exclude empty spaces from being converter to values
                        INVERTER[variable_name] = float(value)
            
            # check if the 'EfficEuroV' (european efficiency) is in the line.
            # If so, stores the european efficiencies in separated variables, each corresponding to a v_eff value
            if 'EfficEuroV' in line:
                key, EfficEuroV_list = line.split('=', 1)  # Split into key and value
                EfficEuroV_list = EfficEuroV_list.strip().split(',')
                for i, value in enumerate(EfficEuroV_list, start=1):
                    variable_name = 'EfficEuroV'+ str(i)
                    if value.strip():      # exclude empty spaces from being converter to values
                        INVERTER[variable_name] = float(value)
            
            # check if the 'ProfilPIOV' is in the line, means that following lines will contain the efficiency values for a voltage
            # line with 'ProfilPIOV1'contain efficiency for lower voltage, 'ProfilPIOV2' for mid voltage and 'ProfilPIOV3' for very high voltage
            elif 'ProfilPIOV' in line:
                V_eff_section = V_eff_section + 1
            
            # if the 'End of TCubicProfile' is in the line and we were in a section where efficiency was outlined, it means the section is over and that voltage points can be resumed in a single variable
            # create the empt y array, sized according to the number of voltage provided
            elif 'End of TCubicProfile' in line and V_eff_section >0:
                INVERTER['p_in_VNomEff' + str(V_eff_section)] = np.zeros(shape = INVERTER.get('NPtsEff'))
                INVERTER['p_out_VNomEff' + str(V_eff_section)] = np.zeros(shape = INVERTER.get('NPtsEff'))
                
                for n in np.arange(1, INVERTER['NPtsEff']+1):
                    key = 'Point_' + str(n)
                    INVERTER['p_in_VNomEff' + str(V_eff_section)][n-1] = INVERTER.get(key)[0]/1000   # converts W in kW
                    INVERTER['p_out_VNomEff' + str(V_eff_section)][n-1] = INVERTER.get(key)[1]/1000  # converts W in kW
            
            # in all other case, just store the variable in the INVERTER dict
            elif '=' in line:  # Look for the presence of '='
                key, value = line.split('=', 1)  # Split into key and value
                value_cleaned = value.strip()  # Clean whitespace
                var_value = _infer_and_convert(value_cleaned)  # Get converted value
                INVERTER[key.strip()] = var_value  # Store only the converted value
    
    # calculate values from Inverter
    
    # remove all unused variables from INVERTER dict
    INVERTER = {key: value for key, value in INVERTER.items() if key in OND_parameters}
    
    # Remap keys to PVlib naming conventions if requested
    if pvlib_remap:
        INVERTER = {pvlib_map.get(k, k): v for k, v in INVERTER.items()}
    
    return INVERTER    


def _infer_and_convert(value: str):
    """
    Infers the data type of a given raw string and converts it accordingly.
    
    It sequentially checks if the string represents a list (comma-separated), 
    a boolean, an integer, or a float. If no specific type is detected, 
    it returns the original string.
    
    Args:
        value (str): The raw string extracted from the text file.
        
    Returns:
        Any: The converted value. Can be a list, bool, int, float, or str.
    """
    
    # 1. List detection (comma-separated sequence)
    # Uses recursion to evaluate and convert each element individually
    if ',' in value:
        value_list = [_infer_and_convert(v.strip()) for v in value.split(',')]
        return value_list
    
    # 2. Boolean detection
    # Handles both True/False and Yes\no common conventions
    lower_value = value.lower()
    if lower_value in {'true', 'false', 'yes', 'no'}:
        return lower_value in {'true', 'yes'}
    
    # 3. Integer detection
    # Validates positive and negative whole numbers
    if value.isdigit() or (value.startswith('-') and value[1:].isdigit()):
        return int(value)
    
    # 4. Float detection
    # Tries to cast to float. Ensures it contains a decimal point or 
    # scientific notation ('e') to avoid overlapping with integers.
    try:
        float_value = float(value)
        if '.' in value or 'e' in lower_value:
            return float_value
    except ValueError:
        pass
    
    # 5. Default fallback
    # If all previous checks fail, return the original string
    return value   

    