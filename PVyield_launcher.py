# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 11:04:49 2024
@author: Antonio Vallone (-T O N I-)
"""

from PVyield import PVyield
from PVyield_classes import Coordinates, Tmy, Power, Mms, Ground, Module, Inverter, MV_Trafo, Field_Curtailment, Loss, MV_Connection, PV_Field, HV_Trafo, HV_Connection, Project_Curtailment, PV_Project
import warnings
from PVyield_classes import custom_formatwarning, WarningLog
from PVyield_report import excel_generator

# Replace the default warning formatter with the custom one
warnings.formatwarning = custom_formatwarning

# Start collecting the warnings, so that they can also be listed in the Excel report.
# This must be run BEFORE the classes below are created, because several warnings are already raised
# while the input classes are being built (e.g. when a PAN or OND file replaces a missing input)
WarningLog.start()

# 0 - PREALLOCATION
excel_generate = True

# ==============================
#        PROJECT CLASSES 
# ==============================
""" Location and TMY sourse are common to all PV felds in PV project """

# 1 - PROJECT LOCATION
coordinates = Coordinates(latitude = 41.902782,     # [deg] Latitude, as decimal degrees (e.g. 20.564) - positive is N emisphere, negative is S emisphere [-90, 90]
                          longitude = 12.496366,     # [deg] Longitude, as decimal (e.g. 15.739) - positive is E emisphere, negative is W emisphere. range [-180, 180]
                          )

# 2 - PROJECT TMY
tmy = Tmy(TMY_source = 'PVGIS',       # [-] TMY source
          api_key = '',               # [-] Only for subscription based services like SolarGIS_API, Solcast_API and SolarAnywhere_API. Default to ''.
          TMY_file_name = '', # [-] Used only for option 'SolarGIS_tmy'. Save the TMY file in the subfolder named TMY inside the same folder where PVsim is located. Default to ''.
          )
"""
# TMY source can either be (NOT case sensitive):
    #'PVGIS'        : automatically download tmy from PVgis v5.3 website for free, using project coordinates
    #'SolarGIS TMY'      : requires to own tmy from solargis. Specify tmy path in variable below
    #'Solcast API'       : requires user to own SolCast API subscription and API key
    #'Meteonorm API'     : requires user to own Meteonorm API subscription and API key
    #'SolarGIS API'      : requires user to own SolarGIS API subscription and API key - TO BE IMPLEMENTED
    #'SolarAnywhere API' : requires user to own Solar Anywhere API subscription and API key - TO BE IMPLEMENTED
"""


# ==============================
#     PV FIELDS
# ==============================
""" Declare one class per type per PV field in the project """

# 3 - PV FIELDS DC AND INVERTER POWER
# Declare one istance per PV field
power_f1 = Power(p_dc_peak = 9.89,         # [MWp] field peak DC power; rounded down to match string count
                 p_ac_nom = 7.5,          # [MVA] field desired nominal inverter power. will be rounded up to match inverter nominal power. to apply curtailment, use Curtailment module
                 )

# 4 - MODULE MOUNTING STRUCTURE
# Declare one istance per field
mms_f1 = Mms (mount_type = 'tracker',    # [-] 'Single-axis tracker' or 'Fixed tilt'. NOT case sensitive
              mod_configuration = '2P',  # [-] '1P', '2P', '3P' ... (Portrait) or '1L', '2L', '3L'... (Landscape). NOT case sensitive
              modules_per_string = 28,   # [-] Number of modules per string in the array.
              pitch = 8.7,               # [m] Row-to-row pitch.
              tilt_angle = 55,           # [deg] MMS tilt/range angle input.
              ground_gap = 0.5,          # [m] Minimum allowed module-to-ground distance (distance between ground and lower border of PVtable). For SAT, is calculated at Max_angle.
              axis_azimuth = -9.3,       # [deg] Compass direction along which the axis of structure lies. Measured clockwise from north (e.g., 90=East, 180=South). Accept also negative values, that will be converter to positive. If None, inferred from mount type (180 for SAT, 90 for FIX). Default to None.
              backtrack = True,          # [-] SAT parameter - backtracking activation. Default to True.)
              )
"""
About tilt_angle:
if SAT, tilt_angle variable is the maximum rotation angle, in decimal degrees, of tracker plane from its
horizontal position (horizontal if tilt_angle = 0). If a float is provided, it represents the maximum rotation angle,
and the minimum rotation angle is assumed to be the opposite of the maximum angle.

If a tuple of (min_angle, max_angle) is provided, it represents both the minimum and maximum rotation angles.
A rotation to 'max_angle' is a counter-clockwise rotation about the y-axis of the tracker coordinate system.
For example, for a tracker with ‘axis_azimuth’ oriented to the south, a rotation to 'max_angle' is towards the west,
and a rotation toward ‘min_angle’ is in the opposite direction, toward the east.
Hence a max_angle of 180 degrees (equivalent to angle = (-180, 180)) allows the tracker to achieve its full rotation capability.

if FIX, tilt_angle variable is the  surface tilt angles in decimal degrees. The tilt angle is defined as degrees from horizontal
(e.g. surface facing up = 0, surface facing horizon = 90)
"""

# 5 - GROUND ALBEDO AND SLOPE 
# Declare one istance per PV field
ground_f1 = Ground(albedo = 0.18,           # [-] Default soil surface albedo. Can be overwritten by imported data (e.g. if using Solcast API).
                   slope_tilt = 0,          # [%] Ground slope in percent.
                   slope_azimuth = 180,     # [deg] Direction of max ground slope gradient. Is the direction where a ball would roll if left on the slope:
                   )                            # 0 = north,  90 = East,  180 = South,  270 = West 
"""
For example, if a site has a slope of 10% in the east direction, the data will be:
slope tilt = 10,   slope_azimuth = 90
"""


# 6 - PV MODULES
# Declare one istance per PV field
module_f1 = Module(PAN_import = True,            # [bool] True if you have a PAN file to import, otherwise set False and then fill manually the data
                   
                   # Automatic MODULE data import from PAN file
                   PAN_file_name = 'Astronergy - CHSM66RN(DG)F-BH-625.PAN',  # [str] name of the PAN file. Save the PAN file in the subfolder named PAN inside the same folder where PVsim is located
                   
                   
                   # Manual MODULE data import - required only if PAN_import = False
                   # The data below are from Trina NEG 690Wp 
                   manufacturer = 'User defined - Trina',  # [-] PV module manufacturer.
                   model = 'NEG 690Wp',      # [-] PV module model.
                   Pnom = 690,               # [Wp] Module DC power at STC.
                   height = 2.384,           # [m] PV module height.
                   width = 1.303,            # [m] PV module width.
                   cells_in_series = 66,     # [-] Number of cells connected in series. Assume full-cut cells equivalent (for half-cut, divide actual number by 2).
                   cells_in_parallel = 2,    # [-] Number of cells in parallel. Typically 2 for half-cut and 1 for full-cut.
                   n_diodes = 3,             # [-] Number of bypass diodes, typically 3.
                   I_L_ref = 18.25,          # [A] short circuit current / Light-generated current (photocurrent) at STC.
                   Voc = 47.9,               # [V] Open circuit voltage at STC.
                   imp_ref = 17.23,          # [A] MPP current at STC.
                   vmp_ref = 40.1,           # [V] MPP voltage at STC.
                   alpha_sc = 6.61,          # [mAmps/Celsius] Short-circuit current temp coefficient. On PAN is 'muIsc' (e.g. 6.57 for Trina NEG 685).
                   gamma_ref = 1.023,        # [-] Diode ideality factor. 1.02 is standard PVsyst model.
                   mu_gamma = -0.0002,       # [1/K] Temperature coefficient for the diode ideality factor.
                   R_sh_ref = 800,           # [ohm] Shunt resistance at STC. If None, inferred from other inputs.
                   R_sh_0 = 80000,           # [ohm] Shunt resistance at zero irradiance. If None, inferred.
                   R_sh_exp = 5.5,           # [-] Single diode exponential parameter. Default a 5.5.
                   R_s = 0.167,              # [ohm] Series resistance at STC. If None, inferred.
                   alpha_absorption = 0.9,   # [-] Module absorptivity (1 - reflectivity).
                   bifaciality = 0.809,      # [-] Bifaciality factor
                   theta_ref = [0, 30, 50, 60, 70, 75, 80, 85, 90],     # [degrees] Array of angles at which the IAM is known.
                   iam_ref = [1, 1, 0.999, 0.994, 0.969, 0.928, 0.829, 0.585, 0],    # [-] IAM values for each angle in iam_theta.
                   V_bypass_diode = 0.7,          # [V] Bypass diode voltage drop. Default a 0.7 (per moduli in silicio).
                   substring_type = 'half cut',   # [-] Either 'half cut' (halfcut cell) or 'full cut' (fullcut cells). NOT case sensitive
                   )


# 7 - INVERTER
# Declare one istance per PV field
inverter_f1 = Inverter(OND_import = True,        # [-] True will import class attribtes from OND file; with False user must define parameters manually.
                       
                       # Automatic INVERTER data import from OND file
                       OND_file_name = 'Huawei - SUN2000-330KTL-H1.OND',    # [-] OND file name. OND file must be saved in 'OND' folder.

                       # Manual INVERTER data import - the variables name matches the one used in OND, dont change them
                       # data below are from Huawei SUN2000-330KTL-H1 string inverter
                       manufacturer = 'User defined - Huawei',  # [-] PV module manufacturer.
                       model = 'SUN2000-330KTL-H1',             # [-] PV module model.
                       Pnom_ac = 300,                           # [kW] nominal AC output power for inverter
                       P_max_ac = 330,                          # [kW] max AC output power for inverter
                       v_min_mppt = 500,                        # [V] min MPPT voltage to allow MPPT to work 
                       v_max_mppt = 1500,                       # [V] max MPPT voltage to allow MPPT to work
                       i_max_mppt = 65,                         # [A] max current per MPPT input
                       n_string_inputs = 28,                    # [-] number of maximum string per inverter
                       n_mppt = 6,                              # [-] number of inverter MPPT
                       v_out = 800,                             # [V] inverter output voltage
                       inv_euro_eff = 98.8,                     # [%] average inverter efficiency (EURO efficiency from datasheet)
                       T_p_nom = 40,                            # [°C] outdoor temperature at which the inverter produces the nominal power - default: 40
                       T_p_max = 30,                            # [°C] outdoor temperature at which the inverter produces the max power - default: 30    
                       )
  

# 8 - LV/MV TRAFO
# Declare one istance per PV field
mv_trafo_f1 = MV_Trafo(v_out = 20,          # [kV] LV/MV trafo output voltage (same as connection voltage if no HV connection is present). Tipically 15, 20 or 30kV
                       p_max = 2.5,         # [MVA] LV/MV transformer max rated power - the simulation will adapt power to field need. defult is 6.4 MVA
                       iron_loss = 0.1,     # [%] iron loss in the LV/MV transformer - default is 0.1%
                       copper_loss = 1.0,   # [%] copper loss in the LV/MV transformer - default is 1.0%
                       )


# 9 - LOSSES
# Declare one istance per PV field (can be also used for multiple PV fields)
loss_f1 = Loss(iam_model = 'GlassAR',          # [-] 'PAN' or 'GlassAR' (NOT case sensitive)
                                                                # 'PAN' will use variables iam_theta and iam_loss in Module object (if available, otherwise will swticth to 'Glass ar')
                                                                # 'GlassAR' will use Martin-ruiz model adjusted for Anti-reflective coated glass
               soiling_loss = 3,               # [%] soiling loss - default is 2%
               rear_shade_factor = 5,          # [%] rear shade factor (fraction of light not hitting the back due to bus and bars). 5% is standard PVsyst assumptions
               thermal_loss_factor = 29,       # [W/m2/K] thermal loss factor for PV module. Default Value is 29, tipical value for most ground mounted applications
               module_quality = -0.3,          # [%] module quality loss (a negative value is a gain, a positive value is a loss) - default is -0.3%
               lid_letid_loss = 0.5,           # [%] LID + LeTID loss - default is 0.5%
               first_year_degrad = 0.3,        # [%] First year degradation - default is 0.3%
               mismatch_loss_front = 0.5,      # [%] array mismatch loss for front cells - default is 0.5%
               mismatch_loss_back = 0.2,       # [%] array mismatch loss for back cells - default is 0.2%
               dc_ohmic_loss_STC = 1.5,        # [%] DC cabling ohmic loss in STC condition - default is 1.5%
               auxiliaries_loss = None,        # [%] Auxiliary loss. If left empty, '', or None, it will be inferred from MMS type
               pv_availability = None,         # [%] PV availability - if None, will infer the value from MMS: 99% for SAT, 99.5 for FT
               ac_lv_cable_len = 500,          # [m] average lenght of the AC LV cable from inverter to trafo. Default is 500m
               ac_lv_cable_sect = 400,         # [mm2] average section of the AC LV cable. Default is 400mm2
               ac_lv_cable_material = 'cu',    # [-] AC LV cable material. Either 'Cu' (copper) or 'Alu' (aluminum). NOT case sensitive. Default is copper
               ac_mv_cable_len = 1000,         # [m] average lenght of the AC MV cable from trafo to collection Cabin. Default is 1000m
               ac_mv_cable_sect = 630,         # [mm2] average section of the AC MV cable from trafo to collection Cabin. DEfault is 630mm2
               ac_mv_cable_material = 'alu',   # [-] AC MV cable material. Either 'Cu' (copper) or 'Alu' (aluminum). NOT case sensitive. Default is alu
               )   

curtailment_f1 = Field_Curtailment(grid_curtail_point = 'Inverter out',  # [-] where apply grid curtail. Available options are: (NOT case sensitive)
                                                                                 # None              : No grid curtailment                                        
                                                                                 # 'Inverter out'    : the grid curtailment is applied at inverter output. AC output power of the inverter is limited to Power.p_ac_inverter for each PV_field. Useful when inverter is curtailed
                                                                                 # 'Field MV cabin'  : the grid curtailment is applied at MV cabin of the field. Useful when PV Field CCI/PPC curtailment is used or PoC is at PV Field collection cabin
                                   grid_curtail_power = 7.5,              # [MVA] grid AC power curtail for point specified)
                                   )

# 10 - MV CONNECTION BETWEEN PV FIELDS
# Declare one istance per PV field
mv_connection_f1 =  MV_Connection(line_len = 2,                 # [km] lenght of the MV connection cable, from collection cabin to HV trafo (if present) or PoC. If None, will assume no MV connection
                                  material = 'alu',             # [str] material of MV connection, either 'Cu' (copper) or 'Alu' (aluminum). Not case sensitive
                                  destination_field = 'grid',   # [str] name of the destination field where the connection is going. Must match an existing PV_Field.name or be 'PoC'/'grid'
                                  auto_sizing= False,           # [bool] if true, number of cable and section will be calculated automatically
                                  line_nr = 1,                  # [-] number of MV connection cable, from  collection cabin to HV trafo (if present) or PoC
                                  line_sect = 630,              # [mm2] section of the MV connection cable, from collection cabin to HV trafo (if present) or PoC
                                  )


# 11 - PV FIELD
# Declare one istance per PV field
field_1 = PV_Field(name = 'PV field #1',                        # [-] name of PV field
                           power = power_f1,                    # object defining DC and AC power for PV field
                           ground = ground_f1,                  # object defining ground orientation of PV field
                           mms =mms_f1,                         # object defining MMS used in PV field          
                           module = module_f1,                  # object defining PV modules used in PV field
                           inverter = inverter_f1,              # object defining inverter used in PV field
                           mv_trafo = mv_trafo_f1,              # object defining MV trafo inside PV field
                           field_curtailment = curtailment_f1,  # object defining field curtailmend for PV field
                           loss = loss_f1,                      # object defining PV field loss
                           mv_connection = mv_connection_f1,    # object defining MV connection outside PV field
                           )

# ====================
# HV CONNECTION
# ====================

# 12 - HV TRAFO
# Declare just one istance common to all pv fields
hv_trafo = HV_Trafo(p_nom = None,          # [MVA] MV/HV transformer rated power. If None, will assume no HV connection
                    iron_loss = None,      # [%] iron loss in the MV/HV transformer - if None, will assign default of 0.4%
                    copper_loss = None,    # [%] copper loss in the MV/HV transformer - if None, will infer data from p_nom
                    )

# 13 - HV CONNECTION
# Declare just one istance common to all pv fields
hv_connection = HV_Connection(line_len = 1)    # [km] lenght of the HV connection cable, from HV trafo (if present) to PoC.
                                               # If None, will assume no HV connection. Loss is 0.1% per km of line

# 14 - PROJECT CURTAILMENT
# Declare just one istance common to all pv fields
curtailment = Project_Curtailment(grid_curtail_point = None,    # [-] where apply grid curtail. Available options are:
                                                              # None                  : No grid curtailment                                        
                                                              # 'Final MV line in'    : the grid curtailment is applied at the entry of the final MV connection line (the one that goes to HV trafo).
                                                              # 'HV trafo in'         : the grid curtailment is applied at the exit of the final MV connection line, i.e. the entry of the HV trafo (if present).
                                                              # 'HV trafo out'        : the grid curtailment is applied at HV trafo output, i.e. the entry of the HV line. Useful when PoC is right after MV/HV trafo
                                                              # 'HV line out'         : the grid curtailment is applied at HV line output. Useful when PoC is at the end of the HV line
                                                              # NOTE: names are NOT case sensitive. 'HV trafo in' is the same physical node as
                                                              #       the exit of the final MV line, and 'HV trafo out' the same as the entry of the HV line.
                                  grid_curtail_power = None,    # [MVA] grid AC power curtail (if any)
                                  )

# ====================
# PV PROJECT
# ====================

pv_project = PV_Project(coordinates = coordinates,          # object defining project location (latitude, longitude)
                        tmy = tmy,                          # object defining TMY source
                        pv_fields = [field_1],              # Tuple/List of PV_Field objects belonging to the project
                        hv_trafo = hv_trafo,                # object defining MV/HV trafo parameters. If None, no HV trafo considered
                        hv_connection =  hv_connection,     # object defining HV connection after MV/HV trafo, If None, no HV connection considered
                        project_curtailment = curtailment,  # object defining power curtailment, if None, no curtailment considered
                        grid_availability = 99.7,           # [%] Grid availability - default is 99.7%
                        )


# ====================
# CHECKS DATA
# ====================
# Every input is validated here. Errors raise an exception and stop the run, while recoverable issues
# raise a warning and the simulation carries on with a standard value
pv_project.input_check()

# ====================
# RUN MODEL
# ====================
# PVyield preallocates the RESULTS DataFrame internally, so there is no need to call Results.assign_results here.
# It returns:
#   RESULTS        the loss chain and the Performance Ratio, one column per PV field plus the 'PV project' column
#   energy_hourly  [kWh] the net energy delivered at the Point of Connection, for each of the 8760 hours of the year
RESULTS, energy_hourly = PVyield(pv_project)

# ====================
# STORE RESULTS
# ====================
if excel_generate:
    # the warnings raised during the run are read automatically from WarningLog and listed in a dedicated sheet
    excel_generator(pv_project, RESULTS, energy_hourly)

# Stop collecting the warnings and put the warnings module back as it was before the run.
# This must always be reached, otherwise the hook installed by WarningLog.start() would stay in place and keep
# collecting warnings raised by any other code run later in the same console session (Spyder, Jupyter, ...).
# The messages collected so far are kept, so WarningLog.get_messages() can still be inspected after the run.
WarningLog.stop()