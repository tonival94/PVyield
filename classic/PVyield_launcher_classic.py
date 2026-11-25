# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 11:04:49 2024
@author: Antonio Vallone (-T O N I-)
"""

from PVyield_classic import PVyield_classic
from PVyield_support_functions_classic import input_check_classic
from PVyield_support_functions_classic import post_processing_classic

# 0 - PREALLOCATION
POWER, LOCATION, METEO, MMS, GROUND, MODULE, INVERTER, TRAFO, LOSS, CONNECTION, CURTAILMENT = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
excel_generate = True

# 1 - DC AND AC POWER
POWER['p_dc_peak'] = 9.892      # [MWp] PV project peak power
POWER['p_ac_connection'] = 7.4  # [MVA] granted connection power as per DNO/TSO

# 2 - LOCATION
LOCATION['latitude'] = 21 #45.2618611  # [deg] Latitude, as decimal degrees (e.g. 20.564) - positive is N emisphere, negative is S emisphere
LOCATION['longitude'] = 79 # 0.4738056 # [deg] Longitude, as decimal (e.g. 20.564) - positive is E emisphere, negative is W emisphere


# 3 - METEO
METEO['TMY_source'] = 'SolarGIS TMY' # [str] TMY source 

# Can either be:
#'PVGIS v5.3'        : automatically download tmy from PVgis v5.3 website for free, using project coordinates
#'SolarGIS TMY'      : requires to own tmy from solargis. Specify tmy path in variable below
#'Solcast API'       : requires user to own SolCast API subscription and API key
#'Meteonorm API'     : requires user to own Meteonorm API subscription and API key
#'SolarGIS API'      : requires user to own SolarGIS API subscription and API key - TO BE IMPLEMENTED
#'SolarAnywhere API' : requires user to own Solar Anywhere API subscription and API key - TO BE IMPLEMENTED

METEO['api_key'] = ''     # [str] only for subcription based services like SolarGIS_API, Solcast_API and SolarAnywhere_API
METEO['TMY_file_name'] = ''  # [str] used only for option 'SolarGIS_tmy'. Save the TMY file in the subfolder named TMY inside the same folder where PVsim is located


# 4 - MODULE MOUNTING STRUCTURE
MMS['mount_type'] = 'Single-axis tracker'  # [str] 'Single-axis tracker' or 'Fixed tilt'
MMS['mod_configuration'] = '2P'            # [str] '1P', '2P', '3P' ... (Portrait) or '1L', '2L', '3L'... (Landscape)
MMS['modules_per_string'] = 28             # [-] Number of modules per string in the array
MMS['pitch'] = 8.7                         # [m] pitch 
MMS['tilt_angle'] = 55                     # [deg] Structure tilt/range angle input - see note below
MMS['axis_azimuth'] = 360-9.3              # [deg] the compass direction along which the axis of structure lies. Measured clockwise from north. (e.g. 90 = east, 180=south, ...). If None, will infer the value from MMS type --> 180 for SAT and 90 for FIX. 
MMS['backtrack'] = True                    # [bool] SAT parameter - backtracking is active as default
MMS['ground_gap'] = 0.5                    # [m] minimum allowed module-to-ground distance. Is the distance between ground and lower border of PVtable. For SAT, is calculated at Max_angle


"""
About tilt_angle:
if SAT, tilt_angle variable is the maximum rotation angle, in decimal degrees, of the one-axis tracker from its
horizontal position (horizontal if axis_tilt = 0). If a float is provided, it represents the maximum rotation angle,
and the minimum rotation angle is assumed to be the opposite of the maximum angle.
If a tuple of (min_angle, max_angle) is provided, it represents both the minimum and maximum rotation angles.
A rotation to 'tilt_angle' is a counter-clockwise rotation about the y-axis of the tracker coordinate system.
For example, for a tracker with ‘axis_azimuth’ oriented to the south, a rotation to 'tilt_angle' is towards the west,
and a rotation toward ‘min_angle’ is in the opposite direction, toward the east.
Hence a max_angle of 180 degrees (equivalent to angle = (-180, 180)) allows the tracker to achieve its full rotation capability.

if FIX, tilt_angle variable is the  surface tilt angles in decimal degrees. The tilt angle is defined as degrees from horizontal
(e.g. surface facing up = 0, surface facing horizon = 90)
"""

# 5 - GROUND 
GROUND['albedo'] = 0.18       # [-] default soil surface albedo. Can be overwritten by imported data (e.g. if using Solast API)
GROUND['slope_tilt'] = 0      # [%] ground slope in percent
GROUND['slope_azimuth'] = 180 # [deg] direction of max ground slope gradient. Is the direction where a ball would roll if left on the slope:
                              # 0 = north
                              # 90 = East
                              # 180 = South
                              # 270 = West               
"""
For example, if a site has a slope of 10% in the east direction, the data will be:
slope tilt = 10
slope_azimuth = 180
"""

# 6 - PV MODULES
# Automatic MODULE data import from PAN file
MODULE['PAN_import'] = True  # [bool] True if you have a PAN file to import, otherwise set False and then fill manually the data
MODULE['PAN_file_name'] = 'Astronergy - CHSM66RN(DG)F-BH-625.PAN' # [str] name of the PAN file. Save the PAN file in the subfolder named PAN inside the same folder where PVsim is located
                               
# Manual MODULE data import - the variables name matches the one used in PAN, dont change them
# The data below are from Trina NEG 690Wp
if not MODULE['PAN_import']:
    MODULE['PNom'] = 690          # [Wp] module DC power at STC
    MODULE['Height'] = 2.384      # [m] PV module height
    MODULE['Width'] = 1.303       # [m] PV module width 
    MODULE['NCelS'] = 66          # [-] number of cells connected in series. This function assume full-cut cells. in halfcut module, cells are half the size and so actual number of cell is double respect to full-cut. In this function, alway putt the full-cut equivalent number of cell (in case of half-cut, take actual number and divide by 2)
    MODULE['NCelP'] = 2           # [-] number of cells in parallel . Tipically =2 in half cut cells and 1 in full cut cells
    MODULE['NDiode'] = 3          # [-] number of Bypass diodes, tipically 3
    MODULE['Isc'] = 18.25         # [A] light-generated current (or photocurrent) at reference conditions, equivalent to short circuit current
    MODULE['Voc'] = 47.9          # [V] open circuit voltage
    MODULE['Imp'] = 17.23         # [A] MPP current at STC condition
    MODULE['Vmp'] = 40.1          # [V] MPP voltage at STC condition
    MODULE['muISC'] = 6.61        # [mAmps/Celsius] short-circuit current temperature coefficient of the module. On PAN is reported as 'muIsc'. For Trina NEG 685 is 6.57 A/C
    MODULE['Gamma'] = 1.023       # [-] diode ideality factor. 1.02 is standard pvsyst model            
    MODULE['muGamma'] = -0.0002   # [1/K] temperature coefficient for the diode ideality factor
    MODULE['RShunt'] = 800        # [ohm] shunt resistance at reference conditions. If None, will infer the data from other inputs
    MODULE['Rp_0'] = 80000        # [ohm] shunt resistance at zero irradiance conditions. If None, will infer the data from other inputs
    MODULE['Rp_Exp'] = 5.5        #[-] Single diode exponential paramenter. Default is 5.5
    MODULE['RSerie'] = 0.167      # [ohm] series resistance at reference conditions. If None, will infer the data from other inputs
    MODULE['Absorb'] = 0.9        # [-] module absorpivity = 1-riflectivity
    MODULE['BifacialityFactor'] = 0.809  # [-] bifaciality factor
    MODULE['iam_theta'] = [0, 30, 50, 60, 70, 75, 80, 85, 90]                 # [degrees] array of angles at which the iam is known 
    MODULE['iam_loss'] = [1, 1, 0.999, 0.994, 0.969, 0.928, 0.829, 0.585, 0]  # [unitless] iam values for each angle in theta_ref 
    MODULE['VRevDiode'] = 0.7     # [V] Bypass diode voltage drop. Default value is 0.7 for Silicium PV modules
    MODULE['SubModuleLayout'] = 'TwinHalfCells'   # [str] either 'half cut' or 'full cit'

# 7 - INVERTER
# Automatic INVERTER data import from OND file
INVERTER['OND_import'] = True  # [bool] True if you have a PAN file to import, otherwise set False and then fill manually the data
INVERTER['OND_file_name'] = 'Huawei - SUN2000-330KTL-H1.OND'   # [str] name of the PAN file. Save the PAN file in the subfolder named PAN inside the same folder where PVsim is located

# Manual INVERTER data import - the variables name matches the one used in OND, dont change them
# data below are from Huawei SUN2000-330KTL-H1 string inverter
if not INVERTER['OND_import']:
    INVERTER['PNomConv'] = 300        # [kW] nominal AC output power of inverter
    INVERTER['PMaxOUT'] = 330         # [kW] max AC output power of inverter
    INVERTER['VMppMin'] = 500         # [V] min MPPT voltage to allow MPPT to work
    INVERTER['VMPPMax'] = 1500        # [V] max MPPT voltage to allow MPPT to work in 
    INVERTER['i_max_mppt'] = 65       # [A] max current per MPPT 
    INVERTER['NbInputs'] = 28         # [-] number of maximum string per inverter
    INVERTER['NbMPPT'] = 6            # [-] number of inverter MPPT
    INVERTER['VOutConv'] = 800        # [V] inverter output voltage 
    INVERTER['inv_euro_eff'] = 98.8   # [%] Average inverter efficiency (EURO efficiency from datasheet)
    INVERTER['TPnom'] = 40            # [°C] outdoor temperature at which the inverter produces the nominal power - default: 40
    INVERTER['TPMax'] = 30            # [°C] outdoor temperature at which the inverter produces the max power - default: 30

# 8 - TRAFO
TRAFO['p_trafo_max'] = 2              # [MVA] max trafo power

# 9 - LOSSES
LOSS['iam_model'] = 'GlassAR' # [str] 'PAN' or 'GlassAR'.
                              # 'PAN' will use variables MODULE['iam_theta'] and MODULE['iam_loss'] or PAN data if provided.
                              # 'GlassAR' will use Martin-ruiz model adjusted for Anti-reflective coated glass
LOSS['soiling_loss'] = 3            # [%] soiling loss - default is 3%
LOSS['rear_shade_factor'] = 5       # [%] rear shade factor (fraction of light not hitting the back due to bus and bars). 5% is standard PVsyst assumptions
LOSS['thermal_loss_factor'] = 29    # [W/m2K] thermal loss factor for PV module. Default Value is 29, tipical value for most ground mounted applications
LOSS['module_quality'] = -0.2       # [%] module quality loss (a negative value is a gain, a positive value is a loss) - default is -0.3%
                                    # NOTE: the standard values assumed by input_check are aligned to the class based version of PVyield
LOSS['lid_loss'] = 0.3              # [%] LID + LeTID loss - default is 0.5%
LOSS['first_year_degrad'] = 0.5     # [%] First year degradation - default is 0.3%
LOSS['mismatch_loss_front'] = 0.6   # [%] array mismatch loss for front cells - default is 0.5%
LOSS['mismatch_loss_back'] = 0.2   # [%] array mismatch loss for back cells - default is 0.2%
LOSS['dc_ohmic_loss_STC'] = 1.0     # [%] DC cabling ohmic loss in STC condition - default is 1.5%
LOSS['trafo_iron_loss'] = 0.1       # [%] iron loss in the LV/MV trasformer - default is 0.1%
LOSS['trafo_copper_loss'] = 1.0     # [%] copper loss in the LV/MV trasformer - default is 1.0%
LOSS['auxiliaries_loss'] = None     # [%] Auxiliary loss. If left empty, '', or None, it will be inferred from MMS type
LOSS['pv_availability'] = None      # [%] PV availability - if None, will infer the value from MMS: 99% for SAT, 99.5 for FT
LOSS['grid_availability'] = 99.7    # [%] Grid availability - default is 99.7%
LOSS['ac_lv_cable_len'] = 250       # [m] average lenght of the AC LV cable from inverter to trafo
LOSS['ac_lv_cable_sect'] = 400      # [mm2] average section of the AC LV cable
LOSS['ac_mv_cable_len'] = 500       # [m] average lenght of the AC MV cable from trafo to collection Cabin
LOSS['ac_mv_cable_sect'] = 240      # [mm2] average section of the AC MV cable from trafo to collection Cabin

# 10 - MV & HV CONNECTION
CONNECTION['mv_v_connection'] = 20            # [kV] trafo connection voltage
CONNECTION['mv_connection_cable_len'] = 0.5   # [km] lenght of the MV connection cable, starting from collection cabin
CONNECTION['mv_connection_cable_sect'] = 185  # [mm2] section of the MV connection cable, starting from collection cabin
CONNECTION['mv_connection_cable_nr'] = 2      # [-] number of MV connection cable, starting from collection cabin

CONNECTION['hv_trafo_p_nom'] = None          # [MVA] MV/HV transformer rated power. If None, will assume no HV connection
CONNECTION['hv_trafo_iron_loss'] = None      # [%] iron loss in the MV/HV transformer - if None, will assign default of 0.4%
CONNECTION['hv_trafo_copper_loss'] = None    # [%] copper loss in the MV/HV transformer - if None, will infer data from p_nom

CONNECTION['hv_connection_cable_len'] = 20    # [km] lenght of the HV connection cable, from HV trafo (if present) to PoC.
                                              # If None, will assume no HV connection. Loss is 0.1% per km of line

# 11 - CURTAILMENT
CURTAILMENT['grid_limitation_point'] = 'POC'   # [-] where to apply grid limitation losses. Available options are:
                                                # 'None'         : no grid limitation applied. The inverter output is not curtailed
                                                # 'inverter'     : the grid limitation is applied at inverter level. useful when inverter is curtailed
                                                # 'MV line in'   : the grid curtailment is applied at the entry of the final MV connection line (the one that goes to HV trafo).
                                                # 'HV trafo in'  : the grid curtailment is applied at the exit of the final MV connection line, i.e. the entry of the HV trafo (if present).
                                                # 'HV trafo out' : the grid curtailment is applied at HV trafo output, i.e. the entry of the HV line. Useful when PoC is right after MV/HV trafo
                                                # 'POC'          : grid limitation applied at POC, after HV connection trafo and line
CURTAILMENT['grid_curtail_power'] = POWER['p_ac_connection']  # [MVA] grid AC power curtail (if any)




# 12 - RUN INPUT CHECK
POWER, LOCATION, METEO, MMS, GROUND, MODULE, INVERTER, TRAFO, LOSS, CONNECTION, CURTAILMENT, errorflag, check_err_msg = input_check_classic(POWER, LOCATION, METEO, MMS, GROUND, MODULE, INVERTER, TRAFO, LOSS, CONNECTION, CURTAILMENT)


# 13 - RUN MODEL
if not errorflag:
    RESULTS, energy_hourly, PVsim_err_msg = PVyield_classic(POWER, LOCATION, METEO, MMS, GROUND, MODULE, INVERTER, TRAFO, LOSS, CONNECTION, CURTAILMENT)
else:
    PVsim_err_msg = '\n\n !! CRITICAL ERROR !! Impossible to proceed with simulation. Check error list and fix inputs'

errors = check_err_msg + PVsim_err_msg
print(errors)

if excel_generate and not errorflag:
    post_processing_classic (RESULTS, energy_hourly, POWER, LOCATION, METEO, MMS, GROUND, MODULE, INVERTER, TRAFO, LOSS, CONNECTION, CURTAILMENT, errors)