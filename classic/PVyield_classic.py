# -*- coding: utf-8 -*-
"""
PVyield classic v13.0
@author: Antonio Vallone (-T O N I-)

Single PV field version of PVyield. It keeps the original structure based on dictionaries, without classes,
and has been aligned to the multi-field version: same bug fixes, same loss chain nomenclature, same
definition of the Performance Ratio and the same Excel report.

Libraries requirements:
> PVlib == 0.16.1
> pandas == 3.0.5
> scipy == 1.18.0
> timezonefinder == 8.5.2
> tzdata == 2026c
> openpyxl == 3.1.5

# improvements for next versions:
# 1 - Add various source for horizon import (e.g. solargis, solcast, etc...) - Waiting for PVlib tomrelease functions
# 2 - Add various API functionality to automatically import TMY from their server (missing: SolarGIS, SolarAnywhere)  - Waiting for PVlib to release functions
# 3 - Add functions to read third party tmy (e.g. solcast). Need a TMY to check how data is stored
# 4 - Add Snow losses calculation (PVlib already includes functions to do that) - snow data can be imported from Solcast API
# 5 - Better horizon loss. Integrate current method in the alternative one
# 6 - Meteonorm TMY import function added, but unable to retrieve data with test API and test location. Need proper API key to be tested
# 7 - Move Grid loss calculation at PoC. Currently this is calculated at inverter output
# 8 - include low irradiation correction in Single diode model. This will fix irradiation losses
# 9 - calculate partial shading loss accurately by performing electric calculation of the series/parallel substrings in the module
# 10 - Add function to automatically calculate cable section and number (from MACRO)
# 11 - Add HV line and HV trafo losses
"""

import os
import math
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import csv
from datetime import datetime
from timezonefinder import TimezoneFinder
import pytz
# import matplotlib.pyplot as plt

from pvlib import location
from pvlib import atmosphere
from pvlib import irradiance
from pvlib import iotools
from pvlib import bifacial
from pvlib import pvsystem
from pvlib import iam
from pvlib import shading
from pvlib import temperature
from pvlib import solarposition
from pvlib import tools
from pvlib import tracking
from pvlib import transformer

from PVyield_support_functions_classic import get_irradiance_full_components
from PVyield_support_functions_classic import inv_Pout_from_Pin
from PVyield_support_functions_classic import ac_ohm_loss
from PVyield_support_functions_classic import vf_row_ground_2d_weights

current_dir = os.path.dirname(os.path.abspath(__file__)) # get current file directory

# Smallest value allowed for a single row view factor. It only acts when the modules lie flat, where the true
# view factor is zero and the corresponding irradiance component is zero as well, turning the ratio used in
# the near shading section into 0/0


# number of structure tilt angles the ground geometry is pre-computed on, then interpolated at every timestep
ground_tilt_steps = 15

def PVyield_classic (POWER, LOCATION, METEO, MMS, GROUND, MODULE, INVERTER, TRAFO, LOSS, CONNECTION, CURTAILMENT):
    
    
    
    # =================================================================================================
    # 1 - UNPACK INPUT VARIABLE
    # =================================================================================================
    
    # 1.0 - HIDDEN VARIABLES AND CONSTANT
    freq = '1h'               # timestep of the simulation. Currently we support only 1h, unless a TMY with higer timestep is imported manually
    horizon_strip_limit = 6.5  # band of the skydome in [deg] where the horizon diffuse component comes from. 6.5° is the standard according to Perex model
    PVsim_err_msg = ''        # inizialize error message to None
    iam_a_r = 0.105           # [-] IAM parameter for martin-ruiz model. Empirical angular loss coefficient. Standard values are 0.08 to 0.25. 0.105 seems to produce best results, in line with Bank advisors. The IAM model is VERY sensible respect to this parameter
    I_o_ref = 0.018E-9        # [A] The dark or diode reverse saturation current at reference conditions, in amperes. 0.018 nA is standard PVsyst value
    EgRef = 1.121             # [eV] The energy bandgap at reference temperature in units of eV. 1.121 eV for crystalline silicon. EgRef must be >0
    irrad_ref = 1000          # [W/m^2] STC irradiance
    temp_ref = 25.0           # [°C] STC temperature
    numerical_zero = 1E-6
    hv_line_loss_per_km = 0.1  # [%] loss of the HV connection line, per km of line
    substrings_per_module = 3 # [-] Number of substring per module
                                  # in Halfcut module, is the substring in half module
                                  # in fullcut module, is the substring in the whole module

    # 1.1 - POWER UNPACK
    p_dc_peak = POWER['p_dc_peak']              # [MWp] PV project peak power
    p_ac_connection = POWER['p_ac_connection']  # [MVA] granted connection power as per DNO
    
    # 1.2 - LOCATION UNPACK    
    latitude = LOCATION['latitude']    # [deg] Latitude, as decimal degrees (e.g. 20.564) - positive is N emisphere, negative is S emisphere
    longitude = LOCATION['longitude']  # [deg] Longitude, as decimal (e.g. 20.564) - positive is E emisphere, negative is W emisphere

    # 1.3 - TMY UNPACK
    TMY_source = METEO['TMY_source'] # [str] TMY source. Can either be:
            #'PVGIS v5.3'        : automatically download tmy from PVgis website for free, using project coordinates
            #'Solcast API'       : requires user to own SolCast API subscription and API key
            #'SolarGIS TMY'      : requires to own tmy from solargis. Specify tmy path in variable below
            #'SolarGIS API'      : requires user to own SolarGIS API subscription and API key - TO BE IMPLEMENTED
            #'SolarAnywhere_API' : requires user to own Solar Anywhere API subscription and API key - TO BE IMPLEMENTED
    api_key = METEO['api_key']   # only for subcription based services like SolarGIS API, Solcast API and SolarAnywhere_API
    TMY_file_name = METEO['TMY_file_name'] # used only for option 'SolarGIS tmy'. Save the TMY file in the subfolder named TMY inside the same folder where PVsim is located

    # 1.4 - MODULE MOUNTING STRUCTURE
    mount_type = MMS['mount_type']                  # [str] 'Single-axis tracker' or 'Fixed tilt'
    mod_configuration = MMS['mod_configuration']    # [str] '1P', '2P', '3P' ... (Portrait) or '1L', '2L', '3L'... (Landscape)
    modules_per_string = MMS['modules_per_string']  # [-] number of modules per string in the array
    pitch = MMS['pitch']                            # [m] pitch
    tilt_angle = MMS['tilt_angle']                  # [deg] Structure tilt/range angle input - seenote below
    axis_azimuth = MMS['axis_azimuth']              # [deg] the compass direction along which the axis of structure lies. Measured east of north. (e.g. 180=south). [degrees] tipically 180 for Single-axis tracker and 90 for Fixed tilt 
    backtrack = MMS['backtrack']                    # [bool] Single-axis tracker parameter - backtracking is active as default
    ground_gap = MMS['ground_gap']                  # [m] minimum allowed module-to-ground distance. Is the distance between ground and lower border of PVtable. For Single-axis tracker, is calculated at Max_angle

    """
    About tilt_angle:
    if Single-axis tracker, tilt_angle variable is the maximum rotation angle, in decimal degrees, of the one-axis tracker from its
    horizontal position (horizontal if axis_tilt = 0). If a float is provided, it represents the maximum rotation angle,
    and the minimum rotation angle is assumed to be the opposite of the maximum angle.
    If a tuple of (min_angle, max_angle) is provided, it represents both the minimum and maximum rotation angles.
    A rotation to 'tilt_angle' is a counter-clockwise rotation about the y-axis of the tracker coordinate system.
    For example, for a tracker with ‘axis_azimuth’ oriented to the south, a rotation to 'tilt_angle' is towards the west,
    and a rotation toward ‘min_angle’ is in the opposite direction, toward the east.
    Hence a max_angle of 180 degrees (equivalent to angle = (-180, 180)) allows the tracker to achieve its full rotation capability.

    if Fixed tilt, tilt_angle variable is the  surface tilt angles in decimal degrees. The tilt angle is defined as degrees from horizontal
    (e.g. surface facing up = 0, surface facing horizon = 90)
    """

    # 1.5 - GROUND UNPACK
    albedo = GROUND['albedo']                # [-] default soil surface albedo. Can be overwritten by imported data (e.g. if using Solast API)
    slope_tilt = GROUND['slope_tilt']        # [%] ground slope in percent
    slope_azimuth = GROUND['slope_azimuth']  # [deg] direction of max ground slope gradient. Is the direction where a ball would roll if left on the slope:0 = north; 90 = East; 180 = South; 270 = West
                    
    """
    For example, if a site has a slope of 10% in the east direction, the data will be:
    slope_azimuth = 180
    slope tilt = 10
    """
    
    # 1.6 - MODULE UNPACK
    # Variables from PAN are remapped to match PVlib nomenclature
    Pnom = MODULE['PNom']                  # [Wp] module DC power at STC
    mod_height = MODULE['Height']           # [m] PV module height
    mod_width = MODULE['Width']             # [m] PV module width 
    cells_in_series = MODULE['NCelS']       # [-] number of cells connected in series. This function assume full-cut cells. in halfcut module, cells are half the size and so actual number of cell is double respect to full-cut. In this function, alway putt the full-cut equivalent number of cell (in case of half-cut, take actual number and divide by 2)
    cells_in_parall = MODULE['NCelP']       # [-] number of cells in parallel . Tipically =2 in half cut cells and 1 in full cut cells
    n_diode =  MODULE['NDiode']             # [-] number of Bypass diodes, tipically 3
    I_L_ref = MODULE['Isc']                 # [A] light-generated current (or photocurrent) at reference conditions, equivalent to short circuit current
    Voc = MODULE['Voc']                     # [V] ppen circuit voltage - USED ONLY DURING INPUT CHECKTO VERIFY MODULE/INVERTER COUPLING
    imp_ref = MODULE['Imp']                 # [A] MPP current at STC condition
    vmp_ref = MODULE['Vmp']                 # [V] The MPP voltage at STC condition
    alpha_sc = MODULE['muISC']              # [mAmps/Celsius] short-circuit current temperature coefficient of the module in units of A/C. On PAN is reported as 'muIsc'. For Trina NEG 685 is 6.57 A/C
    gamma_ref = MODULE['Gamma']             # [-] diode ideality factor. 1.02 is standard pvsyst model            
    mu_gamma = MODULE['muGamma']            # [1/K] temperature coefficient for the diode ideality factor. For Trina NEG 685 is -0.0002 1/K
    R_sh_ref = MODULE['RShunt']             # [ohm] shunt resistance at reference conditions. 800 is stadard PVsyst assumption
    R_sh_0 = MODULE['Rp_0']                 # [ohm] shunt resistance at zero irradiance conditions, in ohms. 80000 is stadard PVsyst assumption
    R_sh_exp = MODULE['Rp_Exp']             # [-] The exponent in the equation for shunt resistance, unitless. Defaults to 5.5
    R_s = MODULE['RSerie']                  # [ohm] series resistance at reference conditions, in ohms.
    alpha_absorption = MODULE['Absorb']     # [-] module absorpivity = 1-riflectivity
    bifaciality = MODULE['BifacialityFactor']  # [-] bifaciality factor
    substring_type = str.lower(MODULE['SubModuleLayout']) # [str] either 'slTwinHalfCells' (halfcut cell) or 'slInLength' (fullcut cell)
    V_bypass_diode = abs(MODULE['VRevDiode'])    # [V] Bypass diode voltage. Default value is 0.7V for Silicium PV modules
    if LOSS['iam_model'] == 'PAN':
        theta_ref = MODULE['iam_theta']     # [deg] array of angles at which the iam is known
        iam_ref = MODULE['iam_loss']        #  [-] iam values for each angle in theta_ref
            
    

    # 1.7 - INVERTER UNPACK
    # Variables from OND are remapped to match PVlib nomenclature
    p_ac_inv_nom = INVERTER['PNomConv']    # [kW] nominal AC output power for inverter
    p_ac_inv_max = INVERTER['PMaxOUT']     # [kW] max AC output power for inverter
    v_min_mppt = INVERTER['VMppMin']       # [V] min MPPT voltage to allow MPPT to work 
    v_max_mppt = INVERTER['VMPPMax']       # [V] max MPPT voltage to allow MPPT to work
    i_max_mppt = INVERTER['i_max_mppt']    # [A] max current per MPPT
    n_inv_inputs = INVERTER['NbInputs']    # [-] number of maximum string per inverter
    n_mppt = INVERTER['NbMPPT']            # [-] number of inverter MPPT
    v_out_inv = INVERTER['VOutConv']       # [V] inverter output voltage
    T_inv_p_nom = INVERTER['TPNom']        # [°C] outdoor temperature at which the inverter produces the nominal power - default: 40
    T_inv_p_max = INVERTER['TPMax']        # [°C] outdoor temperature at which the inverter produces the max power - default: 30
    
    if not INVERTER['OND_import']:
        inv_euro_eff = INVERTER['inv_euro_eff'] # [%] average inverter efficiency (EURO efficiency from datasheet)
        
    else:
        v_eff1 = INVERTER['v_eff1']             # [V] inverter output voltage at which the efficiency is known for p_in_v_eff1 and p_out_v_eff1
        v_eff2 = INVERTER['v_eff2']             # [V] inverter output voltage at which the efficiency is known for p_in_v_eff2 and p_out_v_eff2
        v_eff3 = INVERTER['v_eff3']             # [V] inverter output voltage at which the efficiency is known for p_in_v_eff3 and p_out_v_eff3
        p_in_v_eff1 = INVERTER['p_in_v_eff1']   # [kWp] inverter input DC power at which the inverter efficiency is known for Voltage V_eff1
        p_in_v_eff2 = INVERTER['p_in_v_eff2']   # [kWp] inverter input DC power at which the inverter efficiency is known for Voltage V_eff2
        p_in_v_eff3 = INVERTER['p_in_v_eff3']   # [kWp] inverter input DC power at which the inverter efficiency is known for Voltage V_eff3
        p_out_v_eff1 = INVERTER['p_out_v_eff1'] # [kWac] inverter output AC power resulting from P_in_V_eff1 and V_eff1
        p_out_v_eff2 = INVERTER['p_out_v_eff2'] # [kWac] inverter output AC power resulting from P_in_V_eff2 and V_eff2
        p_out_v_eff3 = INVERTER['p_out_v_eff3'] # [kWac] inverter output AC power resulting from P_in_V_eff3 and V_eff3
        T_eff_1 = INVERTER['TPLim1']            # [°C] outdoor temperature that derates  inverter power to p_out_t1
        T_eff_2 = INVERTER['TPLimAbs']          # [°C] outdoor temperature that derates  inverter power to p_out_t2
        p_out_t1 = INVERTER['PLim1']           # [kW] inverter derated max output power output when outdoor temperature is equal to T_eff_1
        p_out_t2 = INVERTER['PLimAbs']         # [kW] inverter derated max output power output when outdoor temperature is equal to T_eff_2
        

    # 1.8 - TRAFO UNPACK
    p_trafo_max = TRAFO['p_trafo_max']          # [MVA] max trafo power

    # 1.9 - LOSSES UNPACK
    iam_model = LOSS['iam_model']  # [str] 'PAN' or 'GlassAR'.
                                   # 'PAN' will use variables MODULE['iam_theta'] and MODULE['iam_loss'] or PAN data if provided
                                   # 'GlassAR' will use Martin-ruiz model adjusted for Anti-reflective coated glass
    
    soiling_loss = LOSS['soiling_loss']                 # [%] soiling loss - default is 2%
    rear_shade_factor = LOSS['rear_shade_factor']       # [%] rear shade factor (fraction of light not hitting the back due to bus and bars). =5% is standard PVsyst assumptions
    thermal_loss_factor = LOSS['thermal_loss_factor']   # [W/m2K] thermal loss factor for PV module. Default Value is 29, tipical value for most ground mounted applications
    module_quality = LOSS['module_quality']             # [%] module quality loss (a negative value is a gain, a positive value is a loss) - default is -0.2
    lid_loss = LOSS['lid_loss']                         # [%] LID loss - default is 0.5
    first_year_degrad = LOSS['first_year_degrad']       # [%] first year degradation
    mismatch_loss_front = LOSS['mismatch_loss_front']   # [%] array mismatch loss for front cells - default is 1
    mismatch_loss_back = LOSS['mismatch_loss_back']     # [%] array mismatch loss for back cells - default is 0.5
    dc_ohmic_loss_STC = LOSS['dc_ohmic_loss_STC']       # [%] DC cabling ohmic loss in STC condition - default is 1.5
    aux_loss = LOSS['auxiliaries_loss']                 # [%] auxiliary losses - default value is 0.5% for Single-axis tracker and o.3$ for Fixed tilt
    trafo_iron_loss = LOSS['trafo_iron_loss']           # [%] iron loss in the LV/MV trasformer - default is 0.1
    trafo_copper_loss = LOSS['trafo_copper_loss']       # [%] copper loss in the LV/MV trasformer - default is 1
    pv_availability = LOSS['pv_availability']           # [%] PV availability - if None, will infer the value from MMS: 99% for SAT, 99.5 for FT
    grid_availability = LOSS['grid_availability']       # [%] Grid availability - default is 99.7%   
    ac_lv_cable_len = LOSS['ac_lv_cable_len']           # [m] average lenght of the AC LV cable from inverter to trafo
    ac_lv_cable_sect = LOSS['ac_lv_cable_sect']         # [mm2] average section of the AC LV cable [mm2]
    ac_mv_cable_len = LOSS['ac_mv_cable_len']           # [m] average lenght of the AC MV cable from trafos to collection parallel cabin
    ac_mv_cable_sect = LOSS['ac_mv_cable_sect']         # [mm2] average section of the AC MV cable from trafos to collection parallel cabin

    # 1.10 - CONNECTION UNPACK    
    grid_limitation_point = CURTAILMENT['grid_limitation_point'] # [-] where to apply grid limitation losses. Available options are:
                                                                    # 'None'      : no grid limitation applied. The inverter output is not curtailed
                                                                    # '@inverter' : the grid limitation is applied at inverter level. useful when inverter is curtailed
                                                                    # '@POC'      : grid limitation applied at POC. Useful in case a CCI/PPC is used
    grid_curtail_power = CURTAILMENT['grid_curtail_power']      # [MVA] AC power the plant is curtailed to, at the point above
    
    v_connection = CONNECTION['mv_v_connection']                     # [kV] trafo connection voltage
    connection_cable_len = CONNECTION['mv_connection_cable_len']     # [km] lenght of the MV connection cable, from collection parallel cabin onwards
    connection_cable_sect = CONNECTION['mv_connection_cable_sect']   # [mm2] section of the MV connection cable, from collection parallel cabin onwards
    connection_cable_nr = CONNECTION['mv_connection_cable_nr']       # [-] number of MV connection cable, from collection parallel cabin onwards
    
    hv_trafo_p_nom = CONNECTION['hv_trafo_p_nom']                    # [MVA] MV/HV transformer rated power. If None, no MV/HV transformer is present
    hv_trafo_iron_loss = CONNECTION['hv_trafo_iron_loss']            # [%] iron loss of the MV/HV transformer
    hv_trafo_copper_loss = CONNECTION['hv_trafo_copper_loss']        # [%] copper loss of the MV/HV transformer
    hv_connection_cable_len = CONNECTION['hv_connection_cable_len']  # [km] lenght of the HV connection line, from the MV/HV transformer to the PoC
    
    
    
    # =================================================================================================
    # 2 - PREPROCESSING
    # =================================================================================================
    
    # input conversions 
    p_dc_peak = p_dc_peak*1000*1000             # PV project peak power in Wp
    p_ac_connection = p_ac_connection*1000*1000 # granted connection power as per DNO in VA
    
    # the curtailment power is given in MVA like every other power of the launcher, while the whole AC chain
    # works in VA: without this conversion the plant would be clipped to a few VA and the yield would collapse
    if grid_curtail_power is not None:
        grid_curtail_power = grid_curtail_power*1000*1000   # grid curtailment power in VA
    alpha_sc = alpha_sc/1000                    # converts from mAmps/celsius to Amps/celsius
    n_mod = int(mod_configuration[0])           # number of modules in MMS PV table
    
    rear_shade_factor = rear_shade_factor/100   # convert from % to unitless 
    slope_tilt = math.degrees(math.atan(slope_tilt/100)) # convert slope from percentage to degrees
    stc_eff = Pnom/1000/mod_height/mod_width   # STC efficiency of module
    p_ac_inv_nom = p_ac_inv_nom*1000            # inverter AC nominal power in VA
    p_ac_inv_max = p_ac_inv_max*1000            # inverter AC nominal power in VA
    p_trafo_max = p_trafo_max*1000*1000         # max trafo power in VA
    v_connection = v_connection * 1000          # connnection voltage in V
    connection_cable_len = connection_cable_len *1000 # connection cable lenght in m
    is_tracker = 'tracker' in mount_type
    is_portrait = (mod_configuration[1] == 'P')
    is_halfcut = 'half' in substring_type or 'twin' in substring_type
    
   
    if INVERTER['OND_import']:
        p_in_v_eff1 = p_in_v_eff1*1000    # converts from kW to W
        p_in_v_eff2 = p_in_v_eff2*1000    # converts from kW to W
        p_in_v_eff3 = p_in_v_eff3*1000    # converts from kW to W
        p_out_v_eff1 = p_out_v_eff1*1000  # converts from kW to W
        p_out_v_eff2 = p_out_v_eff2*1000  # converts from kW to W
        p_out_v_eff3 = p_out_v_eff3*1000  # converts from kW to W
        p_out_t1 = p_out_t1*1000        # converts from kW to W
        p_out_t2 = p_out_t2*1000        # converts from kW to W
    
    
    
    # =================================================================================================
    # 3 - CREATES RESULTS DATAFRAME
    # =================================================================================================
    
    # The Series holds both numbers and the text of the spacer rows, so it must be of object dtype from the start.
    # It is created with its first row because pandas infers the dtype of an EMPTY Series from the first value added
    # (int64 here), and every later float or text value would then trigger an incompatible-dtype upcast
    RESULTS = pd.Series({'GHI': 0}, dtype=object)
    
    # NOTE: the row names and the block structure are aligned to the multi-field version of PVyield, so that
    # the two loss chains can be compared side by side.
    #
    # The chain is split into three blocks, separated by spacer rows:
    #   BLOCK 1 - the PV plant alone, from GHI down to the exit of the MV cabin
    #   BLOCK 2 - adds the PV availability and the MV connection line leaving the cabin
    #   BLOCK 3 - adds the MV/HV transformer, the HV connection line and the grid availability
    #
    # NOTE: the spacers must be five DIFFERENT strings, otherwise pandas collapses them into a single row and
    # the blocks lose their separation. A different number of trailing spaces makes each of them unique while
    # looking identical once printed
    spacer_1 = '-'*30
    spacer_2 = '-'*30 + ' '
    spacer_3 = '-'*30 + '  '
    spacer_4 = '-'*30 + '   '
    spacer_5 = '-'*30 + '    '
    spacers = [spacer_1, spacer_2, spacer_3, spacer_4, spacer_5]
    
    RESULTS['GHI'] = 0                # will be allocated later with results from calculation
    RESULTS['POA gain'] = 0           # will be allocated later with results from calculation
    
    # ---- BLOCK 1 - the PV plant alone ----
    RESULTS['Horizon loss'] = 0       # will be allocated later with results from calculation
    RESULTS['Near shading loss'] = 0  # will be allocated later with results from calculation
    RESULTS['IAM loss'] = 0           # will be allocated later with results from calculation
    RESULTS['Soiling loss'] = -soiling_loss/100
    RESULTS['POA gain back'] = 0      # will be allocated later with results from calculation
    RESULTS['Temp + Irr loss'] = 0   # will be allocated later with results from calculation
    RESULTS['Partial shading loss'] = 0  # will be allocated later with results from calculation
    RESULTS['Module quality gain'] = -module_quality/100
    RESULTS['LID + LeTID loss'] = -lid_loss/100
    RESULTS['Mismatch front loss'] = -mismatch_loss_front/100
    RESULTS['Mismatch back loss'] = -mismatch_loss_back/100
    RESULTS['DC ohmic wiring loss'] = 0        # will be allocated later with results from calculation
    RESULTS['Inverter MPPT V/I loss'] = 0      # will be allocated later with results from calculation
    RESULTS['Inverter efficiency loss'] = 0    # will be allocated later with results from calculation
    RESULTS['Inverter clipping loss'] = 0      # will be allocated later with results from calculation
    RESULTS['Curtailment @ inverter'] = 0     # will be allocated later with results from calculation
    RESULTS['LV AC cable loss'] = 0           # will be allocated later with results from calculation
    RESULTS['LV/MV trafo loss'] = 0           # will be allocated later with results from calculation
    RESULTS['Auxiliaries loss'] = -aux_loss/100
    RESULTS['MV AC cable loss'] = 0           # will be allocated later with results from calculation
    RESULTS['First year degradation'] = -first_year_degrad/100
    RESULTS['Curtailment @ MV line in'] = 0   # will be allocated later with results from calculation
    RESULTS[spacer_1] = ''                    # spacer
    RESULTS['PR PV only'] = 0                 # will be allocated later with results from calculation
    RESULTS['Yield PV only'] = 0              # will be allocated later with results from calculation
    
    # ---- BLOCK 2 - adds the PV availability and the MV connection line ----
    RESULTS[spacer_2] = ''                    # spacer
    RESULTS['PV availability'] = pv_availability/100-1
    RESULTS['MV connection loss'] = 0         # will be allocated later with results from calculation
    RESULTS[spacer_3] = ''                    # spacer
    RESULTS['PR PV + MV connection + PV avail.'] = 0     # will be allocated later with results from calculation
    RESULTS['Yield PV + MV connection + PV avail.'] = 0  # will be allocated later with results from calculation
    
    # ---- BLOCK 3 - adds the HV parts and the grid availability ----
    RESULTS[spacer_4] = ''                    # spacer
    RESULTS['Curtailment @ HV trafo in'] = 0  # will be allocated later with results from calculation
    RESULTS['HV trafo loss'] = 0              # will be allocated later with results from calculation
    RESULTS['Curtailment @ HV trafo out'] = 0 # will be allocated later with results from calculation
    RESULTS['HV connection loss'] = 0         # will be allocated later with results from calculation
    RESULTS['Curtailment @ PoC'] = 0          # will be allocated later with results from calculation
    RESULTS['Grid availability'] = grid_availability/100-1
    RESULTS[spacer_5] = ''                    # spacer
    RESULTS['PR PV + MV&HV connection + grid avail.'] = 0     # will be allocated later with results from calculation
    RESULTS['Yield PV + MV&HV connection + grid avail.'] = 0  # will be allocated later with results from calculation
    
    
    loss_rows_block_1 = ['Horizon loss', 'Near shading loss', 'IAM loss', 'Soiling loss', 'POA gain back',
                         'Temp + Irr loss', 'Partial shading loss', 'Module quality gain', 'LID + LeTID loss',
                         'Mismatch front loss', 'Mismatch back loss', 'DC ohmic wiring loss',
                         'Inverter MPPT V/I loss', 'Inverter efficiency loss', 'Inverter clipping loss', 'Curtailment @ inverter',
                         'LV AC cable loss', 'LV/MV trafo loss', 'Auxiliaries loss', 'MV AC cable loss',
                         'First year degradation', 'Curtailment @ MV line in']
    
    loss_rows_block_2 = ['PV availability', 'MV connection loss']
    
    loss_rows_block_3 = ['Curtailment @ HV trafo in', 'HV trafo loss', 'Curtailment @ HV trafo out',
                         'HV connection loss', 'Curtailment @ PoC', 'Grid availability']
    
    loss_rows_all = loss_rows_block_1 + loss_rows_block_2 + loss_rows_block_3
    
    pr_rows = ['PR PV only',
               'PR PV + MV connection + PV avail.',
               'PR PV + MV&HV connection + grid avail.']
    
    yield_rows = ['Yield PV only',
                  'Yield PV + MV connection + PV avail.',
                  'Yield PV + MV&HV connection + grid avail.']
    
    # set units for the results above.
    results_units = pd.Series(index = RESULTS.index, dtype=object)
    results_units.loc[:] = '%'
    results_units.loc['GHI'] = 'kWh/m2'
    results_units.loc[yield_rows] = 'kWh/kWp'
    results_units.loc[spacers] = ''
    
    
    
    # =================================================================================================
    # 4 - GENRATE TIMESTEPS FOR SIMULATION
    # =================================================================================================
    
    # identify current timezone from coordinates
    timezone_name = TimezoneFinder().timezone_at(lat=latitude,
                                                 lng=longitude,
                                                 )
    
    # Get the UTC offset
    timezone = pytz.timezone(timezone_name)
    UTC_offset = timezone.utcoffset(datetime.now()).total_seconds() / 3600  # Offset in hours

    
    local_time = pd.date_range(start ='1990-01-01 00:00', 
                               end ='1990-12-31 23:00',
                               freq = freq,
                               tz=timezone,
                                )
    
    
    
    # =================================================================================================
    # 5 - SET LOCATION
    # =================================================================================================
    try:
        altitude = location.lookup_altitude(latitude = latitude,
                                            longitude = longitude,
                                            )
    except:
        PVsim_err_msg += 'WARNING! Not able to retrieve altitude for given location. The simulation assumed an altitude of 20m above sea level/n/n'
        altitude = 20
    
    site_location = location.Location(latitude = latitude,
                                      longitude = longitude,
                                      tz = timezone,
                                      name = None,
                                      altitude = altitude,
                                      )
    
    
    
    # =================================================================================================
    # 6 - IMPORT RADIATION DATA
    # =================================================================================================
    
    # 6.1 - IMPORT TMY DATA FROM SELECTED SOURCE
    
    try:
        if TMY_source == 'PVGIS v5.3':
            imported_tmy, *_ = iotools.get_pvgis_tmy(latitude = latitude,
                                                     longitude = longitude, 
                                                     outputformat = 'csv', 
                                                     usehorizon = False, 
                                                     userhorizon = None, 
                                                     startyear = 2005, 
                                                     endyear = 2015, 
                                                     map_variables = True, 
                                                     url = 'https://re.jrc.ec.europa.eu/api/', 
                                                     timeout=30,
                                                     )
            
        elif TMY_source == 'Solcast API':
        # To test the API with no paid subscription, register to solcast and then use coordinates:
        # latitude =  -33.856784
        # longitude = 151.215297
            imported_tmy, *_ = iotools.get_solcast_tmy(latitude = latitude,
                                                       longitude = longitude, 
                                                       api_key = api_key,
                                                       map_variables = True,
                                                       output_parameters = ['ghi','dni','dhi','air_temp','surface_pressure','relative_humidity','wind_speed_10m','albedo'],
                                                       # Other output can be specified as specified in https://docs.solcast.com.au/#output-parameter-definitions
                                                       )
            
        elif TMY_source == 'Meteonorm API':
        # To test the API with no paid subscription, register to solcast and then use coordinates:
        # latitude =  50
        # longitude = 10
        # API_key = 'demo0000-0000-0000-0000-000000000000'
            imported_tmy, *_ = iotools.get_meteonorm_tmy(latitude = latitude,
                                                         longitude = longitude,
                                                         api_key = api_key,
                                                         parameters = 'all',
                                                         time_step = '1h',
                                                         horizon = 'flat',
                                                         map_variables = True,
                                                         # Other output can be specified as specified in https://docs.meteonorm.com/api/climate#climate-tmy
                                                         )
            
            # the PVlib map_variables of Meteonorm does not cover the air pressure, so it keeps its original Meteonorm name and its original unit hPa.
            # It is renamed and converted to Pa here. NOTE: The parameter name is looked up rather than hardcoded, because it has changed across
            # Meteonorm API versions. If no pressure column is found, the standard value is used further below
            for meteonorm_pressure_name in ['surface_air_pressure', 'air_pressure', 'pressure', 'atmospheric_pressure']:
                if meteonorm_pressure_name in imported_tmy:
                    imported_tmy = imported_tmy.rename(columns={meteonorm_pressure_name: 'pressure'})
                    imported_tmy['pressure'] = imported_tmy['pressure'] * 100
                    break
            
        # The online databases above deliver the TMY in UTC. The simulation runs in local time,
        # therefore the data are rolled by the UTC offset of the site.
        # NOTE: the roll is applied column by column so that the dtype of each column is preserved
        if TMY_source in ['PVGIS v5.3', 'Solcast API', 'Meteonorm API']:
            imported_tmy = pd.DataFrame(np.roll(imported_tmy, shift=int(round(UTC_offset)), axis=0), columns=imported_tmy.columns)
            imported_tmy.index = local_time
            
    # if unable to import data from online databases above, the data are st to 0
    except:
        PVsim_err_msg += 'ERROR! Not able to get radiation data from online database. Not able to compute simulation\n\n'
        imported_tmy = pd.DataFrame(index = local_time)
        imported_tmy['ghi'] = 0
        imported_tmy['dhi'] = 0
        imported_tmy['dni'] = 0
        imported_tmy['temp_air'] = 25
        imported_tmy['pressure'] = 101325
        imported_tmy['relative_humidity'] = 0
        imported_tmy['wind_speed'] = 0
            
 
    
    if TMY_source == 'SolarGIS TMY':
        TMY_path = os.path.join(current_dir, 'TMY', TMY_file_name)
        # finds the rows where data TMY data starts
        with open(TMY_path) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=';')
            line_count = 1
            for row in csv_reader:
                if row[0] == '#Data:':
                    break
                else:
                    line_count = line_count+1
    
    
        # read TMY specified in user path
        imported_tmy = pd.read_csv(TMY_path,
                                   sep=';',
                                   header=line_count,
                                   )
        # SolarGIS TMY are already in local time, no need to shift data to match our timesteps
        imported_tmy.index = local_time
    
        # remap variables to match PVlib nomenclature
        imported_tmy = imported_tmy.rename(columns={'GHI':  'ghi',
                                                    'DIF':  'dhi',
                                                    'DNI':  'dni',
                                                    'TEMP': 'temp_air',
                                                    'AP':   'pressure',
                                                    'RH' :  'relative_humidity',
                                                    'WS':   'wind_speed'}
                                                     )
        
        # SolarGIS TMY reports the air pressure (AP) in hPa, while every PVlib function expects it in Pa.
        if 'pressure' in imported_tmy:
            imported_tmy['pressure'] = imported_tmy['pressure'] * 100
    
    # 6.2 UNPACK DATA INTO LOCAL VARIABLES
    ghi = imported_tmy['ghi']
    dhi = imported_tmy['dhi']
    dni = imported_tmy['dni']
    temp_air = imported_tmy['temp_air']
    wind_speed = imported_tmy['wind_speed']
    
    # check if optional variable have been imported (not all TMY import methods support these variables)
    if 'albedo' in imported_tmy:
        albedo = imported_tmy['albedo']   # overwrites the albedo values with the one in the TMY
    
    # A sanity check is also added: every source is converted to Pa where needed, so a pressure of a few
    # hundred here means a provider changed its unit or a new source was added without its own conversion
    if 'pressure' in imported_tmy:
        pressure = imported_tmy['pressure']
        if np.nanmedian(pressure) < 50000:
            PVsim_err_msg += 'WARNING! The air pressure imported from the TMY is far below the 101325 Pa expected at sea level. Data are probably in hPa rather than Pa. The simulation used the standard pressure instead\n\n'
            pressure = pd.Series(101325.0, index=local_time, name="air pressure")
    else:
        pressure = pd.Series(101325.0, index=local_time, name="air pressure")  # if no pressure present in TMY, average air pressure is used instead
    
    if 'relative_humidity' in imported_tmy:
        relative_humidity = imported_tmy['relative_humidity'] # at the monment not used
    
    # Stores GHI into results Dataframe
    RESULTS['GHI'] = np.nansum(ghi)/1000



    # =================================================================================================
    # 7 - COMPUTE TIME SHIFT
    # =================================================================================================
    """
    In the PVlib convention, the time label always refers to the beginning of the record, and concerns the accumulation up to the next record
    If times = 7:00, means the radiation data are the sum of radiation from 7:00 to 8:00, and solar positions are computed at 7:00
    Now, with imported data files, there is no general rule.
    Each record usually holds a time label, which possibly doesn't correspond to the PVlib definition.
    Moreover, this time label may represent data measured (averaged)  over a time period which doesn't correspond to the true Legal Time "round" hours
    We need to align the imported TMY with our timesteps. This is done by trying to shift the imported data and see what shift gives a best match with clear sky irradiation
    """
    
    # for TMY source SolarGIS TMY and Solcast API the timeshift is knows 
    if TMY_source in ['SolarGIS TMY', 'Solcast API']:
        time_shift = pd.Timedelta(minutes = 30)
    
    else:
        # Create a vector with all timeshift that will be tried, from -30 min to + 30 min
        err = pd.Series(index = np.arange(-30,31,1))
        
        # cycles through 60 min timeshift (+/- 30 min) and compute the Mean Square Error for each tentative shift
        for time_shift in np.arange(-30,30,1):
            time_delta = pd.Timedelta(minutes=time_shift)
            
            clear_sky = site_location.get_clearsky(times = local_time + time_delta,
                                                   model='ineichen',
                                                   )
            # assign to clear_sky the same index as ghi to allow for numpy calculation with ghi
            # NOTE: clear_sky is computed with the timeshift applied, we simply reassign the original index so that we can preform operations with ghi variable
            clear_sky.index = local_time
            
            # computes Mean Square Error and stores it
            MSE = np.power(clear_sky['ghi']-imported_tmy['ghi'],2).mean()
            err.loc[time_shift] = MSE 
        
        
        time_shift = err.idxmin()                       # find the time_shift associated to the minimum MSE
        time_shift = pd.Timedelta(minutes=time_shift)   # convert the timeshift in a pandas timedelta
    
    # compute the new shifted time steps that will be used for solar geometry calculations
    local_time_calc = local_time + time_shift          
    
    
    
    # =================================================================================================
    # 8 - CALCULATE SOLAR GEOMETRY
    # =================================================================================================
    
    # compute solar positions (azimuth and zenith)
    solar_position = solarposition.get_solarposition(time = local_time_calc,
                                                     latitude = latitude,
                                                     longitude = longitude,
                                                     altitude = altitude,
                                                     pressure = pressure,
                                                     method = 'nrel_numpy',
                                                     temperature = temp_air,
                                                     )
    
    solar_position.index = local_time # change to local time for better readability
    
    
    # compute extraterrestrial normal irradiance
    dni_extra = irradiance.get_extra_radiation(datetime_or_doy = local_time_calc,
                                               solar_constant = 1361,
                                               method ='spencer',
                                               epoch_year = datetime.now().year,
                                               )
    
    dni_extra.index = local_time # change UTC time for better readability
    
    #computes relative airmass
    airmass_relative = atmosphere.get_relative_airmass(zenith = solar_position['apparent_zenith'],
                                                       model='gueymard2003',
                                                       )
    
    
    
    # =================================================================================================    
    # 9 - ALLOCATE PV STRUCTURE DATA AND GET ARRAY ORIENTATION
    # =================================================================================================
    
    # compute basic geometric data
    if is_portrait:
        pvrow_width = mod_height * n_mod #width of PVtable when Portrait, including 2mm bezel
    elif not is_portrait:
        pvrow_width = mod_width * n_mod  #width of PVtable when Landascape, including 2mm bezel
    
    pvrow_height = ground_gap + pvrow_width/2 * np.sin(np.radians(tilt_angle)) # Height of the pvtable, measured at the tablerow midline (equivalent to Torque tube axis for Single-axis tracker), tipically 1.5m for 1V and 2.5m for 2V (m)
    gcr = pvrow_width / pitch
       
    
    # Allocate mount (Single-axis tracker or Fixed tilt) data structure for PVlib functions
    
    # Calculate tracker axis tilt in the global reference frame when on a sloped plane.
    # Axis tilt is the inclination of the tracker rotation axis with respect to horizontal, ranging from 0 degrees (horizontal axis) to 90 degrees (vertical axis)
    axis_tilt = tracking.calc_axis_tilt(slope_azimuth = slope_azimuth,
                                        slope_tilt = slope_tilt,
                                        axis_azimuth = axis_azimuth,
                                        )
    
    # Calculate the angle, relative to horizontal, of the line formed by the intersection
    # between the slope containing the tracker axes and a plane perpendicular to the tracker axes.
    # See details here: https://pvlib-python.readthedocs.io/en/stable/reference/generated/pvlib.shading.shaded_fraction1d.html
    cross_axis_tilt = tracking.calc_cross_axis_tilt(slope_azimuth = slope_azimuth,
                                                    slope_tilt = slope_tilt,
                                                    axis_azimuth = axis_azimuth, 
                                                    axis_tilt = axis_tilt,
                                                    )
    
    if is_tracker:
       max_angle = tilt_angle
       mount = pvsystem.SingleAxisTrackerMount(axis_tilt = axis_tilt,
                                               axis_azimuth = axis_azimuth,
                                               max_angle = max_angle,
                                               backtrack = backtrack,
                                               gcr = gcr,
                                               racking_model = 'open_rack',
                                               module_height = pvrow_height,
                                               cross_axis_tilt = cross_axis_tilt,
                                               )
    
    else:
       surface_tilt = tilt_angle
       surface_azimuth = axis_azimuth+90  # relevant for Fixed tilt only. Surface azimuth angles in decimal degrees. Surface_azimuth must be >=0 and <=360. The Azimuth convention is defined as degrees east of north (e.g. North = 0, South=180 East = 90, West = 270)
       mount = pvsystem.FixedMount(surface_tilt = surface_tilt,
                                     surface_azimuth = surface_azimuth,
                                     racking_model='open_rack',
                                     module_height = pvrow_height,
                                     )
       
    # compute structure orientation data, like surface tilt, azimuth, aoi
    orientation = mount.get_orientation(solar_zenith = solar_position['apparent_zenith'],
                                        solar_azimuth = solar_position['azimuth'],
                                        )
    
    # Compute angle of incidence for front and back
    aoi_front = irradiance.aoi(surface_tilt = orientation ['surface_tilt'],
                               surface_azimuth = orientation['surface_azimuth'],
                               solar_zenith = solar_position['apparent_zenith'],
                               solar_azimuth = solar_position['azimuth'])
    
    aoi_back = 180-aoi_front
    
    
    
    # =================================================================================================
    # 10 - GET POA IRRADIATION FOR A SINGLE SHED USING PEREZ MODEL (SAME AS PVSYST)
    # =================================================================================================
    
    irradPOA_single_shed = get_irradiance_full_components(surface_tilt = orientation['surface_tilt'],
                                                          surface_azimuth = orientation['surface_azimuth'],
                                                          solar_zenith = solar_position['apparent_zenith'],
                                                          solar_azimuth = solar_position['azimuth'],
                                                          aoi = aoi_front,
                                                          ghi = ghi,
                                                          dhi = dhi,
                                                          dni = dni,
                                                          dni_extra = dni_extra,
                                                          airmass = airmass_relative,
                                                          albedo = albedo,
                                                          )
    
    # compute POA gain over front side
    RESULTS['POA gain'] = np.nansum(irradPOA_single_shed['poa_global'])/np.nansum(ghi)-1 
    
    
    
    # =================================================================================================
    # 11 - GET HORIZON DATA & LOSSES
    # =================================================================================================
    
    # 11.1 - GET HORIZON DATA
    # get horizon data from PVGIS. Result is a tuple
    try:
        hor_data = iotools.get_pvgis_horizon(latitude = latitude,
                                            longitude = longitude,
                                            url = 'https://re.jrc.ec.europa.eu/api/',
                                            )

    # if horizon cannot be retireved, assume flat horizon
    except:
        hor_data = [pd.Series(data = [0,0], index = [0,180], name = 'horizon_elevation'),'Not able to retrieve horizon data from PVGIS. The simulation assumed flat horizon']
        PVsim_err_msg += ('WARNING! Could not get horizon elevation data from PVGIS. The simulation assumed no horizon')

    # extract data from tuple and creates Dataframe to store future horizon data
    hor_data = pd.DataFrame(hor_data[0], index = hor_data[0].index)
    
    
    # 11.2 - GET IDEAL SUNRISE SUNSET TIMES    
    #for each day, compute ideal sunrise and sunset time, corresponding to flat horizon
    day_times = pd.date_range(start ='1990-01-01 00:00',   
                              end ='1990-12-31 23:59',
                              freq = '1D',
                              tz = timezone,
                              )
    
    sun_rise_set_year = solarposition.sun_rise_set_transit_spa(times = day_times + pd.Timedelta(days=0.5),   # by specifing the hour as noon (+0.5days respect to original timestep)), the function sunrise_set_transit get the correc day of the sunrise and sunset. ih midnight is passed, the sunrise is calculated on previous day
                                                               longitude = longitude,
                                                               latitude = latitude)
    
    # rename the Dataframe colums for better readability and change index to the one matching the rest of PV lib calculations 
    sun_rise_set_year = sun_rise_set_year.rename(columns={'sunrise':'ideal_sunrise', 'sunset':'ideal_sunset'})
    sun_rise_set_year = sun_rise_set_year.drop('transit', axis = 1)
    sun_rise_set_year.index = day_times

    # 11.3 - GET REAL SUNRISE SUNSET TIMES
    # real sunrise and sunset are computed fromimported horizon data
    for day in day_times:       
        
        precise_time = pd.date_range(start = day,
                                     end = day + pd.Timedelta(hours=23, minutes=59),
                                     freq='1min',
                                     tz = timezone,
                                     )
        
        # get the current hour of precise_time. To avoid issues with DST (Day saving time) we must first decolasize the data, then set minutes and seconds to zero, then relocalize the data
        current_hour = precise_time.tz_localize(None).floor('h').tz_localize(timezone, ambiguous = False)   
        
        # compute precise sun position for every minute. This will be used to get the precise time when the sun cross the horizon
        sun_positions_rise_set = solarposition.get_solarposition(time = precise_time,
                                                                 latitude = latitude,
                                                                 longitude = longitude,
                                                                 altitude = altitude,
                                                                 pressure = pressure.loc[current_hour],
                                                                 method = 'nrel_numpy',
                                                                 temperature = temp_air.loc[current_hour],
                                                                 )
        
        # determine sun apparent elevation corresponsing to the azimuth where horizon profile is provided
        hor_elev = np.interp(sun_positions_rise_set['azimuth'], hor_data.index, hor_data['horizon_elevation'])
        sun_above_horizon = np.array(sun_positions_rise_set['apparent_elevation'] > hor_elev)
    
        # get effective sunrise and sunset time, accounting for horizon
        sun_rise_set_year.loc[day, 'real_sunrise'] = pd.Series(precise_time).where(sun_above_horizon, pd.NaT).dropna().min()  # sunrise time is the earliest time when the sun is above the horizon
        sun_rise_set_year.loc[day, 'real_sunset'] = pd.Series(precise_time).where(sun_above_horizon, pd.NaT).dropna().max()  # sunset time is the last time sun is above the horizon
        
    # 11.4 - GET DNI HORIZON LOSS
    dni_loss_hor = pd.Series(index = local_time)
    
    for hours in local_time:
        hour_start = hours
        hour_end = hour_start + pd.Timedelta(hours=1)
        current_day = hours.normalize()
        
        if hour_end <= sun_rise_set_year.loc[current_day, 'real_sunrise'] or hour_start >= sun_rise_set_year.loc[current_day, 'real_sunset']:      # if timestep ends before sunrise or starts after sunset, the DNI is lost
            dni_loss_hor.loc[hours] = 1 
        elif hour_start >= sun_rise_set_year.loc[current_day, 'real_sunrise'] and hour_end <= sun_rise_set_year.loc[current_day, 'real_sunset']:   # if the timestamp start and finish are within the sunrse and sunset time, no DNI loss at all
            dni_loss_hor.loc[hours] =  0
        elif hour_start < sun_rise_set_year.loc[current_day, 'real_sunrise'] < hour_end:
            dni_loss_hor.loc[hours] = (sun_rise_set_year.loc[current_day, 'real_sunrise'] - hour_start) / pd.Timedelta(hours=1)
        elif hour_start < sun_rise_set_year.loc[current_day, 'real_sunset'] < hour_end:
            dni_loss_hor.loc[hours] = (hour_end - sun_rise_set_year.loc[current_day, 'real_sunset']) / pd.Timedelta(hours=1)
        
    dni_hor = dni * (1-dni_loss_hor) 
      
    
    # 11.5 - DHI HORIZON LOSS
    # Compute the delta across the horizon azimuthal data (it is the 'Base' of the integration shape used later).
    # Note this process will create an array with one value less then hor_elev_angles.
    # Hence it is added a final value to restore removed data and complete the full rotation 0-360°
    hor_data['delta_azimuth'] = np.append( np.diff(hor_data.index), 360 - hor_data.index[-1])  
    
    # compute the area of the projected area of the horizon over the skydome by integrating all the horizon pieces
    # the area of a rectangular sphere sector in steradians from equator is delta_azimuth * (1- cos(elevation_angle))
    # since the area is not rectangular but trapeziodal, we will calculate the area using the right and left elevation, then we will average
    area_horizon_strip_rx = np.sum(hor_data['delta_azimuth']/360*2*np.pi * (1-tools.cosd(hor_data['horizon_elevation'])))
    area_horizon_strip_lx = np.sum(hor_data['delta_azimuth']/360*2*np.pi * (1-tools.cosd(np.roll(hor_data['horizon_elevation'],1))))
    
    area_horizon_strip = (area_horizon_strip_rx+area_horizon_strip_lx)/2
    
    # compute the area of the sky dome in steradians (is the suface of half a sphere)
    area_sky_dome = 2 * np.pi
    
    # compute dhi horizon losses are the ration between the ara of the horizon over the area of the sky dome
    dhi_hor = dhi * (1-area_horizon_strip/area_sky_dome)
    
    # 11.6 - COMPUTE POA HORIZON LOSSES
    # recompute ghi from dhi_hor and dhi_hor
    ghi_hor = dhi_hor + dni_hor * tools.cosd(solar_position['zenith'])
        
    irradPOA_hor = get_irradiance_full_components(surface_tilt = orientation['surface_tilt'],
                                                     surface_azimuth = orientation['surface_azimuth'],
                                                     solar_zenith = solar_position['apparent_zenith'],
                                                     solar_azimuth = solar_position['azimuth'],
                                                     aoi = aoi_front,
                                                     ghi = ghi_hor,
                                                     dhi = dhi_hor,
                                                     dni = dni_hor,
                                                     dni_extra = dni_extra,
                                                     airmass = airmass_relative,
                                                     albedo = albedo,
                                                     )
        
    RESULTS['Horizon loss'] = np.nansum(irradPOA_hor['poa_global'])/np.nansum(irradPOA_single_shed['poa_global'])-1 
    
    
    
    """
    # ALTERNATIVE HORIZON LOSS METHOD -- DO NOT DELETE
    # This method calculate horizon loss on each POA component, insted of calculating it on the horizontal components
    # could be included in get_irradiance_full_components function to easily calculate POA front and POA back irradiance with and without horizon
    
    # preallocate Dataframe to store irradiation data after horizon losses
    irradPOA_hor = pd.DataFrame(index = irradPOA_single_shed.index, columns = irradPOA_single_shed.columns)
    
    # 11.2 - HORIZON LOSS FOR DIRECT & CIRCUMSOLAR COMPONENT
    hor_data_direct = np.interp(solar_position['azimuth'], hor_data.index, hor_data['horizon_elevation'])  # Interpolate the horizon elevation data to the solar azimuth and keep as a numpy array.
    hor_data_direct = pd.Series(hor_data_direct, local_time)                                               # Convert to Pandas Series to allow usage in PVlib functions
    
    # compute direct and circumsolar irradiation data after horizon loss
    irradPOA_hor['poa_direct'] = irradPOA_single_shed['poa_direct'].where(solar_position['apparent_elevation'] > hor_data_direct, 0) # NOTE: in Pandas.where, when the condition is True, the original value is kept. If condition is False, the new value is replaced with alternative value
    irradPOA_hor['poa_sky_diffuse_circumsolar'] = irradPOA_single_shed['poa_sky_diffuse_circumsolar'].where(solar_position['apparent_elevation'] > hor_data_direct, 0) # NOTE: in Pandas.whre, when the condition is True, the original value is kept. If condition is False, the new value is replaced with alternative value
    
    irradPOA_hor['poa_direct'] = irradPOA_single_shed['poa_direct'] * (1-dni_loss_hor)
    irradPOA_hor['poa_sky_diffuse_circumsolar'] = irradPOA_single_shed['poa_sky_diffuse_circumsolar'] * (1-dni_loss_hor)
    
    
    # 11.3 - HORIZON LOSS FOR HORIZON DIFFUSE COMPONENTS
    # to compute the horizon loss for horizon diffuse component, we must calculate what percentage of horizon sky dome strip is covered by the horizon profile
    
    # Compute the delta across the horizon azimuthal data (it is the 'Base' of the integration shape used later).
    # Note this process will create an array with one value less then hor_elev_angles.
    # Hence it is added a final value to restore removed data and complete the full rotation 0-360°
    hor_data['delta_azimuth'] = np.append( np.diff(hor_data.index), 360 - hor_data.index[-1]) 
    
    #  define the horizon boundary limit. Tipically it is assumed to between 5 and 15°. Here we assume 10°. This is the sky dome strip where the diffuse horizon irradiation comes from
    # Cap the horizon profile to the horizon boundary limit [0,horizon_strip_limit]
    hor_data['horizon_elevation_hor_strip'] = np.clip(hor_data['horizon_elevation'], 0, horizon_strip_limit)
    
    # compute area of horizon profile for diffuse horizon losses in steradians
    # the area of a rectangular sphere sector from equator is delta_azimuth * (1- cos(elevation_angle)) * R^2. R does not matter since we will consider the ration bwtween two areas
    # since the area is not rectangular but trapeziodal, we will calculate the area using the right and left elevation, then we will average
    area_horizon_strip_rx = np.sum(hor_data['delta_azimuth']/360*2*np.pi * (1-tools.cosd(hor_data['horizon_elevation_hor_strip'])))
    area_horizon_strip_lx = np.sum(hor_data['delta_azimuth']/360*2*np.pi * (1-tools.cosd(np.roll(hor_data['horizon_elevation_hor_strip'],1))))
    
    area_horizon_strip = (area_horizon_strip_rx+area_horizon_strip_lx)/2
    
    # compute the area of the sky dome strip between 0° e horizon_strip_limit (in steradians), where diffuse horizon irradiation comes from
    area_horizon_limit =2*np.pi * (1-tools.cosd(horizon_strip_limit))
    
    # calculate diffuse horizon loss as the ratio between horizon strip and the portion of the sky associated with horizon irradiation
    irradPOA_hor['poa_sky_diffuse_horizon'] = irradPOA_single_shed['poa_sky_diffuse_horizon']*(1-area_horizon_strip/area_horizon_limit)
    
    
    # 11.4 - HORIZON LOSS FOR SKY ISOTROPIC DIFFUSE COMPONENTS
    # to compute the horizon loss for sky isotropic component, we must calculate what percentage of sky dome is covered by the horizon profile
    hor_data['horizon_elevation_sky_strip'] = np.where(hor_data['horizon_elevation']>horizon_strip_limit, hor_data['horizon_elevation'], horizon_strip_limit)
    
    # compute area of horizon profile for diffuse sky losses in steradians
    # the area of a rectangular sphere sector from equator is delta_azimuth * (1- cos(elevation_angle)) * R^2. R does not matter since we will consider the ration bwtween two areas
    # since the area is not rectangular but trapeziodal, we will calculate the area using the right and left elevation, then we will average
    area_sky_strip_rx = np.sum(hor_data['delta_azimuth']/360*2*np.pi * (1-tools.cosd(hor_data['horizon_elevation_sky_strip'])))
    area_sky_strip_lx = np.sum(hor_data['delta_azimuth']/360*2*np.pi * (1-tools.cosd(np.roll(hor_data['horizon_elevation_sky_strip'],1))))
    
    # only the portion of horizon profile exceeding the horizon strip limit is accounted
    area_sky_strip = (area_sky_strip_rx+area_sky_strip_lx)/2 - area_horizon_limit
    
    # calculate the area of the sky dome (in steradians) where sky isotropic irradiation comes from. Equal to 2π (full sky dome = half of a sphere) minus horizon limit strip
    total_sky_area = 2 * np.pi - area_horizon_limit
    
    irradPOA_hor['poa_sky_diffuse_isotropic'] = irradPOA_single_shed['poa_sky_diffuse_isotropic']*(1-area_sky_strip/total_sky_area)
    
    # 11.5 COMPUTE REST OF RADIATION COMPONENTS AND LOSSES
    irradPOA_hor['poa_sky_diffuse'] = irradPOA_hor['poa_sky_diffuse_isotropic'] + irradPOA_hor['poa_sky_diffuse_circumsolar'] + irradPOA_hor['poa_sky_diffuse_horizon']
    irradPOA_hor['poa_ground_diffuse'] = irradPOA_single_shed['poa_ground_diffuse']
    irradPOA_hor['poa_diffuse'] = irradPOA_hor['poa_sky_diffuse'] + irradPOA_hor['poa_ground_diffuse']
    irradPOA_hor['poa_global'] = irradPOA_hor['poa_direct'] + irradPOA_hor['poa_diffuse']
    
    RESULTS['Horizon loss'] = np.nansum(irradPOA_hor['poa_global'])/np.nansum(irradPOA_single_shed['poa_global'])-1 
    """
    
    
    # =================================================================================================
    # 12 - BACK IRRADIANCE
    # =================================================================================================
    
    irradPOA_back = get_irradiance_full_components(surface_tilt = 180 - orientation['surface_tilt'],
                                                    surface_azimuth = (180 + orientation['surface_azimuth']),
                                                    solar_zenith = solar_position['apparent_zenith'],
                                                    solar_azimuth = solar_position['azimuth'],
                                                    aoi = aoi_back,
                                                    ghi = ghi_hor,
                                                    dhi = dhi_hor,
                                                    dni = dni_hor,
                                                    dni_extra = dni_extra,
                                                    airmass = airmass_relative,
                                                    albedo = albedo,
                                                    )
    
    
    
    # =================================================================================================
    # 13 - SHADING FACTORS
    # =================================================================================================
    
    # 13.1 - SHADING FACTOR - DIRECT RADIATION - FRONT AND BACK
    # calculate shaded fraction of module row width due to self shading
    SF_direct_front = shading.shaded_fraction1d(solar_zenith = solar_position['apparent_zenith'],
                                                solar_azimuth = solar_position['azimuth'],
                                                axis_azimuth = axis_azimuth,
                                                shaded_row_rotation = orientation['surface_tilt'],
                                                collector_width = pvrow_width,
                                                pitch = pitch,
                                                axis_tilt = axis_tilt,
                                                surface_to_axis_offset = 0.1,  # distance between torque tube and module surface. For Fixed tilt, is the distance beteen support structure plane and module surface 
                                                cross_axis_slope = cross_axis_tilt,
                                                shading_row_rotation = orientation['surface_tilt'],
                                                )
    
    SF_direct_back = shading.shaded_fraction1d(solar_zenith = solar_position['apparent_zenith'],
                                               solar_azimuth = solar_position['azimuth'],
                                               axis_azimuth = axis_azimuth,
                                               shaded_row_rotation = 180 - orientation['surface_tilt'],
                                               collector_width = pvrow_width,
                                               pitch = pitch,
                                               axis_tilt = axis_tilt,
                                               surface_to_axis_offset = 0.1,  # distance between torque tube and module surface. For Fixed tilt, is the distance beteen support structure plane and module surface 
                                               cross_axis_slope = cross_axis_tilt,
                                               shading_row_rotation = 180-orientation['surface_tilt'],
                                               )
    
    
    # 13.2 - SHADING FACTORS - DIFFUSE SKY RADIATION - FRONT AND BACK
    # since we have POA component, to compute shading loss, we need to know the difference of view factor (module->sky) from the single shed to the multi shed case
    # the ratio among the two will give the loss of POA irradiance due to multi shed case
    
    # calculate view Factor from module surface to sky, for front and back side, for both one row and multiple row case
    VF_PV_sky_front_multi_row = bifacial.utils.vf_row_sky_2d_integ(surface_tilt = orientation['surface_tilt'],
                                                                   gcr = gcr,
                                                                   x0=0,
                                                                   x1=1,
                                                                   )
    
    # calculate view Factor from module surface to sky, for front side, for both one row and multiple row case
    VF_PV_sky_back_multi_row = bifacial.utils.vf_row_sky_2d_integ(surface_tilt = 180 - orientation['surface_tilt'],
                                                                  gcr = gcr,
                                                                  x0=0,
                                                                  x1=1,
                                                                  )
    
    
    VF_PV_sky_front_single_row = (1 + tools.cosd(orientation['surface_tilt']))/2   # never zero for a tilt below 180 deg
    
    # View factor for sky is clipped to zero because formula can give slighly negative values for some angles instead of zero.
    # Numerical zero is used instead of zero because view factor is denominator in future fomrmulas. A value equal to zero would resul in cose error
    if is_tracker:
        VF_PV_sky_back_single_row = np.maximum((1 + tools.cosd(180-orientation['surface_tilt']))/2, numerical_zero)
    else:
        VF_PV_sky_back_single_row = max((1 + tools.cosd(180-orientation['surface_tilt']))/2, numerical_zero)
    
    # remove the negative values PVlib returns at the extremes of the tilt range, which are rounding artefacts
    VF_PV_sky_front_multi_row = np.maximum(VF_PV_sky_front_multi_row, 0)
    VF_PV_sky_back_multi_row = np.maximum(VF_PV_sky_back_multi_row, 0)
    
    
    # 13.3 - SHADING FACTOR - DIFFUSE GROUND RADIATION - FRONT AND BACK
    """
    Since we have POA component, to compute shading loss, we need to know the difference of view factor (module->ground) from the single shed to the multi shed case
    Respect to single shed case, there is also a reduced radiation from sky to ground, since the module will prevent the ground to receive diffuse radiation.
    We will acount this reduction with the View factor Sky-ground for multi shed case VF_gnd_sky_multi_row
    
    Also, differently from single shed case, the direct-like radiation could also be blocked by the module!
    Hence we need to compute the portion of the ground that receives diffuse radiation only (because shaded by the modules) and the portion receiving direct+circumsolar+diffuse
    
    Therefore, we also need to know  also need to know the circumsolar + direct component sum on the horizontal plane
    """
    
    # compute the PV row to consider in integration in order to have a integration step less than 5°
    max_rows = np.ceil(pvrow_height / (pitch * tools.tand(5)))

    # compute the view factor from ground to sky in case of multiple sheds.
    # NOTE: since section 8.6.4 resolves the brightness of the ground point by point, this pitch-averaged value is no longer used by the calculation. 
    # NOTE: since PVlib 0.16 the integral is computed analytically (parameters npoints and vectorize were removed)
    # and the parameter surface_tilt was renamed tracker_rotation. PVlib itself passes the surface tilt here.
    VF_gnd_sky_multi_row = bifacial.utils.vf_ground_sky_2d_integ(tracker_rotation = orientation['surface_tilt'],
                                                                 gcr = gcr,
                                                                 height = pvrow_height,
                                                                 pitch = pitch,
                                                                 max_rows = max_rows,
                                                                 )
    
    # calculate view Factor from module surface to ground, for front and back side, for both one row and multiple row case
    VF_PV_ground_front_multi_row = bifacial.utils.vf_row_ground_2d_integ(surface_tilt = orientation['surface_tilt'],
                                                                         gcr = gcr, 
                                                                         x0=0,
                                                                         x1=1,
                                                                         )
    
    VF_PV_ground_back_multi_row = bifacial.utils.vf_row_ground_2d_integ(surface_tilt = 180 - orientation['surface_tilt'],
                                                                        gcr = gcr, 
                                                                        x0=0,
                                                                        x1=1,
                                                                        )
    
    # View factor for ground is clipped to zero because formula can give slighly negative values for some angles instead of zero.
    # Numerical zero is used instead of zero because view factor is denominator in future fomrmulas. A value equal to zero would resul in cose error
    if is_tracker:
        VF_PV_ground_front_single_row = np.maximum((1 - tools.cosd(orientation['surface_tilt']))/2, numerical_zero)
    else:
        VF_PV_ground_front_single_row = max((1 - tools.cosd(orientation['surface_tilt']))/2, numerical_zero)
    
    # the back side looks at (1+cos(tilt))/2, which is never zero for a tilt below 180 deg, hence no floor needed
    VF_PV_ground_back_single_row  = (1 - tools.cosd(180 - orientation['surface_tilt']))/2
    
    
    # compute the fraction of ground between rows that receives direct radiation, accounting for the shade of the panels.
    # NOTE: as for VF_gnd_sky_multi_row above, this pitch-averaged value is no longer used by the calculation. A copy is kept for future usages
    # Since PVlib 0.16 this private function takes the projected solar angle phi instead of the solar and surface angles,
    # so phi is computed first, the same way PVlib does in bifacial.infinite_sheds
    phi_gnd = tools.atand(bifacial.utils._solar_projection_tangent(solar_zenith = solar_position['apparent_zenith'],
                                                                   solar_azimuth = solar_position['azimuth'],
                                                                   surface_azimuth = orientation['surface_azimuth'],
                                                                   ))
    fraction_gnd_direct = bifacial.utils._unshaded_ground_fraction(tracker_rotation = orientation['surface_tilt'],
                                                                   phi = phi_gnd,
                                                                   gcr = gcr,
                                                                   )
    
    # compute circumsolar and direct component on the horizontal plane
    # to do so, we assumed a flat surface       
    irrad_horizontal = get_irradiance_full_components(surface_tilt = 0,
                                                           surface_azimuth = orientation['surface_azimuth'],
                                                           solar_zenith = solar_position['apparent_zenith'],
                                                           solar_azimuth = solar_position['azimuth'],
                                                           aoi = solar_position['apparent_zenith'],  # for a flat surface, the incidence angle is equal to the sun zenith angle
                                                           ghi= ghi_hor,
                                                           dhi = dhi_hor,
                                                           dni = dni_hor,
                                                           dni_extra = dni_extra,
                                                           airmass = airmass_relative,
                                                           albedo = 0,   # we dont need this since it is already accounted and we only need the ratio between direct+circumsolar and diffuse
                                                           )
    
    direct_fraction  = (irrad_horizontal['poa_direct'] + irrad_horizontal['poa_sky_diffuse_circumsolar'])/irrad_horizontal['poa_global']
    diffuse_fraction = 1 - direct_fraction
    
    
    
    # =================================================================================================
    # 13.1 - GROUND BRIGHTNESS SEEN BY THE MODULES
    # =================================================================================================
    # The ground has very non uniform brightness; in fact some areas (e.g. under the module) are shadowed, while other areas (between rows) are lit by sun.
    # By calculating the average of ground brightness, I would overestimate the ground component.
    # Hence, ground component is calcualted by brightness weighted averaged of each ground section with its brightness. View factor and ground brighness are caluclated for each ground strip and then integrated all together.
    # However, for tracker the view factor an brightness chances at each angle. therefore, we calculate such data for 15 standard angle. All other tracker position are calcualted via interpolation of such 7 standard positions
    # This is the same method used by PVsyst

    # Compute the standard tilt angles on which the geometry is pre-computed. A tracker sweeps its whole
    # range, for a fixed structure a single angle is enough
    if is_tracker:
        standard_angles = np.linspace(0, np.nanmax(orientation['surface_tilt']), ground_tilt_steps)
    else:
        standard_angles = np.array([float(np.nanmax(orientation['surface_tilt']))])
    
    # weights of the view of the ground from each face, for every standard angle at once
    x_ground, weights_ground_front, weights_ground_back = vf_row_ground_2d_weights(standard_angles = standard_angles,
                                                                                   gcr = gcr,
                                                                                   pvrow_height = pvrow_height,
                                                                                   pitch = pitch,
                                                                                   )
    
    # Pre-compute the sky seen from every ground point. It does not depend on the face, only on the structure.
    # PVlib returns the result shaped (ground points, angles), so it is transposed to match the weights above
    vf_sky_ground = np.asarray(bifacial.utils.vf_ground_sky_2d(rotation = standard_angles,
                                                               gcr = gcr,
                                                               x = x_ground,
                                                               pitch = pitch,
                                                               height = pvrow_height,
                                                               max_rows = 10,
                                                               )).T
    
    # interpolate the pre-computed quantities on the actual tilt of every timestep.
    # For a fixed structure the tilt is a constant, for a tracker it is a time series: np.broadcast_to turns
    # both into an array of the right length without duplicating any memory
    tilt_hourly = np.broadcast_to(np.asarray(orientation['surface_tilt'], dtype=float), (len(local_time),))
    tilt_filled = np.nan_to_num(tilt_hourly, nan=0.0)
    
    if len(standard_angles) == 1:
        # fixed structure: the same geometry applies to every timestep
        weights_front_hourly = np.broadcast_to(weights_ground_front[0], (len(tilt_filled), len(x_ground)))
        weights_back_hourly = np.broadcast_to(weights_ground_back[0], (len(tilt_filled), len(x_ground)))
        vf_sky_hourly = np.broadcast_to(vf_sky_ground[0], (len(tilt_filled), len(x_ground)))
    else:
        # tracker: linear interpolation between the two standard angles bracketing the actual tilt.
        # It is done by hand rather than with np.interp because that only works on one dimension at a time,
        # and here a whole ground profile has to be interpolated at every timestep
        lower = np.clip(np.searchsorted(standard_angles, tilt_filled) - 1, 0, len(standard_angles)-2)
        weight_upper = np.clip((tilt_filled - standard_angles[lower])/(standard_angles[lower+1] - standard_angles[lower]), 0, 1)[:, np.newaxis]
        
        weights_front_hourly = weights_ground_front[lower]*(1-weight_upper) + weights_ground_front[lower+1]*weight_upper
        weights_back_hourly = weights_ground_back[lower]*(1-weight_upper) + weights_ground_back[lower+1]*weight_upper
        vf_sky_hourly = vf_sky_ground[lower]*(1-weight_upper) + vf_sky_ground[lower+1]*weight_upper
    
    # solar projection angle on the cross-axis plane, positive when the sun is on the side the front faces.
    # It places the shadow of the rows on the ground at every timestep
    tan_psi = np.asarray(bifacial.utils._solar_projection_tangent(solar_zenith = solar_position['apparent_zenith'],
                                                                  solar_azimuth = solar_position['azimuth'],
                                                                  surface_azimuth = orientation['surface_azimuth'],
                                                                  ), dtype=float)
    tan_psi = np.nan_to_num(tan_psi, nan=0.0, posinf=0.0, neginf=0.0)
    
    # THE SHADOW OF THE ROWS ON THE GROUND, resolved at every timestep.
    # Each edge of the row is projected onto the ground along the direction of the sun: a point at height h
    # lands h*tan_psi further along the ground. With the right edge of the row higher, the front looks to the
    # left, so a sun on the front side casts the shadow towards the right.
    # All the quantities below are expressed as a fraction of the pitch, like x_ground
    cos_tilt, sin_tilt = tools.cosd(tilt_filled), tools.sind(tilt_filled)
    
    # shadow_1: ground point where the RIGHT (higher) edge of the row is projected
    shadow_1 = (gcr/2*pitch*cos_tilt + (pvrow_height + gcr/2*pitch*sin_tilt)*tan_psi)/pitch
    
    # shadow_2: ground point where the LEFT (lower) edge of the row is projected
    shadow_2 = (-gcr/2*pitch*cos_tilt + (pvrow_height - gcr/2*pitch*sin_tilt)*tan_psi)/pitch
    
    # the shadow of one row is the band between those two points, described by its centre and its half width.
    # An extra axis is added so that they broadcast against the ground positions further below
    shadow_centre = ((shadow_1 + shadow_2)/2)[:, np.newaxis]
    shadow_half_width = (np.abs(shadow_1 - shadow_2)/2)[:, np.newaxis]
    
    # lit_hourly: 1 where a ground point receives the beam, 0 where it falls inside the shadow of a row.
    # The rows being infinite and regularly spaced, the shadow pattern repeats every pitch, so the distance
    # from the centre of the shadow is taken modulo one pitch before comparing it with the half width
    lit_hourly = (np.abs(((x_ground[np.newaxis, :] - shadow_centre + 0.5) % 1.0) - 0.5) > shadow_half_width).astype(float)
    
    # brightness of every ground point, then averaged with the weights of the view of each face
    brightness_hourly = (np.asarray(diffuse_fraction, dtype=float)[:, np.newaxis]*vf_sky_hourly
                         + np.asarray(direct_fraction, dtype=float)[:, np.newaxis]*lit_hourly)
    
    ground_brightness_front = pd.Series(np.einsum('ij,ij->i', weights_front_hourly, brightness_hourly), index=local_time)
    ground_brightness_back = pd.Series(np.einsum('ij,ij->i', weights_back_hourly, brightness_hourly), index=local_time)



    # =================================================================================================
    # 14 - NEAR SHADING LOSSES
    # =================================================================================================
    
    # create struct to store irradPOA_near_shade (front+back)
    POA_columns = ['poa_global',
                   'poa_front',
                   'poa_front_direct',
                   'poa_front_diffuse',
                   'poa_front_ground_diffuse',
                   'poa_front_sky_diffuse',
                   'poa_front_sky_diffuse_isotropic',
                   'poa_front_sky_diffuse_circumsolar',
                   'poa_front_sky_diffuse_horizon',
                   'poa_back',
                   'poa_back_direct',
                   'poa_back_diffuse',
                   'poa_back_ground_diffuse',
                   'poa_back_sky_diffuse',
                   'poa_back_sky_diffuse_isotropic',
                   'poa_back_sky_diffuse_circumsolar',
                   'poa_back_sky_diffuse_horizon']
    
    irradPOA_near_shade = pd.DataFrame(columns = POA_columns)
    
    # POA direct radiation after near shading front and back
    irradPOA_near_shade['poa_front_direct'] = irradPOA_hor['poa_direct'] * (1-SF_direct_front)
    irradPOA_near_shade['poa_back_direct'] = irradPOA_back['poa_direct'] * (1-SF_direct_back) * (1-rear_shade_factor)
    
    # POA sky isotropic diffuse radiation after near shading front and back
    irradPOA_near_shade['poa_front_sky_diffuse_isotropic'] = irradPOA_hor['poa_sky_diffuse_isotropic'] / VF_PV_sky_front_single_row * VF_PV_sky_front_multi_row 
    irradPOA_near_shade['poa_back_sky_diffuse_isotropic'] = irradPOA_back['poa_sky_diffuse_isotropic'] / VF_PV_sky_back_single_row  * VF_PV_sky_back_multi_row * (1-rear_shade_factor)
    
    # POA sky circumsolar diffuse radiation after near shading front and back. We assume circumsolar behaves like direct radiation, so shading factor is applied
    irradPOA_near_shade['poa_front_sky_diffuse_circumsolar'] = irradPOA_hor['poa_sky_diffuse_circumsolar'] * (1-SF_direct_front)
    irradPOA_near_shade['poa_back_sky_diffuse_circumsolar'] = irradPOA_back['poa_sky_diffuse_circumsolar'] * (1-SF_direct_back) * (1-rear_shade_factor)
    
    # POA sky horizon diffuse radiation after near shading
    # this is set to zero assuming that modules are not able to see horizon strip. In fact horizon strip bad is 6.5° according to Perez. Such band is likely hidden by next PV structure
    # in previous caluclation, it was assumed that horizon view factor was the same as sky view factor. However, such assumption is completely wrong and produced overestimated horizon irradiance
    irradPOA_near_shade['poa_front_sky_diffuse_horizon'] = 0  # irradPOA_near_shade['poa_front_sky_diffuse_horizon'] = irradPOA_hor['poa_sky_diffuse_horizon'] / VF_PV_sky_front_single_row * VF_PV_sky_front_multi_row
    irradPOA_near_shade['poa_back_sky_diffuse_horizon'] = 0  # irradPOA_near_shade['poa_back_sky_diffuse_horizon'] = irradPOA_back['poa_sky_diffuse_horizon'] / VF_PV_sky_back_single_row * VF_PV_sky_back_multi_row * (1-rear_shade_factor)
    
    # compute summary of diffuse components
    irradPOA_near_shade['poa_front_sky_diffuse'] = irradPOA_near_shade['poa_front_sky_diffuse_isotropic'] + irradPOA_near_shade['poa_front_sky_diffuse_circumsolar'] + irradPOA_near_shade['poa_front_sky_diffuse_horizon']
    irradPOA_near_shade['poa_back_sky_diffuse'] = irradPOA_near_shade['poa_back_sky_diffuse_isotropic'] + irradPOA_near_shade['poa_back_sky_diffuse_circumsolar'] + irradPOA_near_shade['poa_back_sky_diffuse_horizon']
    
    # POA ground diffuse radiation after near shading, front and back.
    # The chain is: (irradiance leaving the open ground) x (ground brightness reduction) x (view factor row -> visible ground)
    #   - dividing by the single row view factor gives back the irradiance leaving an unobstructed ground (ghi*albedo)
    #   - the multi row view factor accounts for how much of that ground the module actually sees
    #   - ground_brightness front and back is an averaged reduction factor that trasforms irradiance leaving unobstructed ground (ghi*albedo) to the actual irradiance leaving the ground under the module
    
    irradPOA_near_shade['poa_front_ground_diffuse'] = irradPOA_hor['poa_ground_diffuse'] / VF_PV_ground_front_single_row * VF_PV_ground_front_multi_row * ground_brightness_front
    irradPOA_near_shade['poa_back_ground_diffuse'] =  irradPOA_back['poa_ground_diffuse'] / VF_PV_ground_back_single_row * VF_PV_ground_back_multi_row * ground_brightness_back * (1-rear_shade_factor)
    
    # PREVIOUS FORMULATION - DO NOT DELETE
    # Previously, instead of using the averaged ground_brightness, it was used (diffuse_fraction * VF_gnd_sky_multi_row + direct_fraction * fraction_gnd_direct)
    #   - diffuse_fraction: fraction of diffuse radiation in POA radiation (=1-direct_fraction)
    #   - VF_gnd_sky_multi_row: View factor from ground to sky
    #   - direct_fraction: fraction of direct radiation in POA radiation (=1-diffuse_fraction)
    #   - fraction_gnd_direct: fraction of area under PV structures that receives direct radiation
    # irradPOA_near_shade['poa_front_ground_diffuse'] = irradPOA_hor['poa_ground_diffuse'] / VF_PV_ground_front_single_row * VF_PV_ground_front_multi_row * (diffuse_fraction * VF_gnd_sky_multi_row + direct_fraction * fraction_gnd_direct)
    # irradPOA_near_shade['poa_back_ground_diffuse'] = irradPOA_back['poa_ground_diffuse'] / VF_PV_ground_back_single_row * VF_PV_ground_back_multi_row * (diffuse_fraction * VF_gnd_sky_multi_row + direct_fraction * fraction_gnd_direct) * (1-rear_shade_factor)

    
    # complete the irradPOA_near_shade variable with remaining components.
    irradPOA_near_shade['poa_front_diffuse'] = irradPOA_near_shade['poa_front_sky_diffuse'] + irradPOA_near_shade['poa_front_ground_diffuse']
    irradPOA_near_shade['poa_back_diffuse'] = irradPOA_near_shade['poa_back_sky_diffuse'] + irradPOA_near_shade['poa_back_ground_diffuse']
    irradPOA_near_shade['poa_front'] = irradPOA_near_shade['poa_front_direct'] + irradPOA_near_shade['poa_front_diffuse']
    irradPOA_near_shade['poa_back'] = irradPOA_near_shade['poa_back_direct'] + irradPOA_near_shade['poa_back_diffuse']
    irradPOA_near_shade['poa_global'] = irradPOA_near_shade['poa_front'] + irradPOA_near_shade['poa_back']
    
    # to respect PVsyst convenction, the shading loss is compute don front only.
    # However, the dataframe irradPOA_near_shade contains already the POA for back side including shade loss
    # Such components will continue to be used and the POA gain on back side will be calculated later
    RESULTS['Near shading loss'] = np.nansum(irradPOA_near_shade['poa_front'])/np.nansum(irradPOA_hor['poa_global'])-1
    
    
    
    # =================================================================================================
    # 15 - IAM LOSS
    # =================================================================================================
    
    # 15.1 - DIRECT COMPONENT IAM
    # of the IAM model is Glass or Glass AR, the IAM angles are generated and its relative IAM values are computed
    if iam_model =='GlassAR':
        theta_ref = np.arange(0, 90, 0.5)      
        iam_ref = iam.martin_ruiz(aoi = theta_ref,
                                  a_r = iam_a_r,
                                  )
        
    # IMA direct compotentn loss is calculated as interpolation of the angle of incidence over the IAM loss function
    iam_direct_front = iam.interp(aoi = aoi_front,
                                  theta_ref = theta_ref,
                                  iam_ref = iam_ref, 
                                  method='cubic', 
                                  normalize=True,
                                  )
        
    iam_direct_back = iam.interp(aoi = aoi_back,
                                 theta_ref = theta_ref,
                                 iam_ref = iam_ref,
                                 method='cubic', 
                                 normalize=True,
                                 )

        
    # 15.2 - DIFFUSE COMPONENT IAM
    #compute IAM loss function for diffuse radiation
    iam_function = interp1d(theta_ref, iam_ref, kind='cubic', bounds_error=False, fill_value=(iam_ref[0], iam_ref[-1]))
    
    
    # compute IAM los for each separate diffuse component
    iam_diffuse_front_sky = iam.marion_integrate(function = iam_function,
                                                 surface_tilt = np.max(orientation['surface_tilt']),
                                                 region = 'sky',
                                                 num=None,
                                                 )
        
    iam_diffuse_front_horizon = iam.marion_integrate(function = iam_function,
                                                     surface_tilt = np.max(orientation['surface_tilt']),
                                                     region = 'horizon',
                                                     num=None,
                                                     )
        
    iam_diffuse_front_ground = iam.marion_integrate(function = iam_function,
                                                    surface_tilt = np.max(orientation['surface_tilt']),
                                                    region = 'ground',
                                                    num=None,
                                                    )
        
    iam_diffuse_back_sky = iam.marion_integrate(function = iam_function,
                                                surface_tilt = 180 - np.max(orientation['surface_tilt']),
                                                region = 'sky',
                                                num=None,
                                                )
        
    iam_diffuse_back_horizon = iam.marion_integrate(function = iam_function,
                                                    surface_tilt = 180 - np.max(orientation['surface_tilt']),
                                                    region = 'horizon',
                                                    num=None,
                                                    )
    
    iam_diffuse_back_ground = iam.marion_integrate(function = iam_function,
                                                   surface_tilt = 180 - np.max(orientation['surface_tilt']),
                                                   region = 'ground',
                                                   num=None,
                                                   )
    
    
    # preallocate a new dictionary for irradPOA_iam, then calculate components after IAM loss
    irradPOA_iam = pd.DataFrame(columns = POA_columns)
    # IAM loss for direct radiation
    irradPOA_iam['poa_front_direct'] = irradPOA_near_shade['poa_front_direct'] * iam_direct_front    
    irradPOA_iam['poa_back_direct'] = irradPOA_near_shade['poa_back_direct'] * iam_direct_back
    
    # IAM for sky diffuse radiation components. For circumsolar, we apply IAM loss of direct component
    irradPOA_iam['poa_front_sky_diffuse_isotropic'] = irradPOA_near_shade['poa_front_sky_diffuse_isotropic'] *iam_diffuse_front_sky
    irradPOA_iam['poa_back_sky_diffuse_isotropic'] = irradPOA_near_shade['poa_back_sky_diffuse_isotropic'] *iam_diffuse_back_sky
    irradPOA_iam['poa_front_sky_diffuse_circumsolar'] = irradPOA_near_shade['poa_front_sky_diffuse_circumsolar'] * iam_direct_front
    irradPOA_iam['poa_back_sky_diffuse_circumsolar'] = irradPOA_near_shade['poa_back_sky_diffuse_circumsolar'] * iam_direct_back
    irradPOA_iam['poa_front_sky_diffuse_horizon'] = irradPOA_near_shade['poa_front_sky_diffuse_horizon'] * iam_diffuse_front_horizon
    irradPOA_iam['poa_back_sky_diffuse_horizon'] = irradPOA_near_shade['poa_back_sky_diffuse_horizon'] * iam_diffuse_back_horizon
    
    # compute summary of sky diffuse components
    irradPOA_iam['poa_front_sky_diffuse'] = irradPOA_iam['poa_front_sky_diffuse_isotropic'] + irradPOA_iam['poa_front_sky_diffuse_circumsolar'] + irradPOA_iam['poa_front_sky_diffuse_horizon']
    irradPOA_iam['poa_back_sky_diffuse'] = irradPOA_iam['poa_back_sky_diffuse_isotropic'] + irradPOA_iam['poa_back_sky_diffuse_circumsolar'] + irradPOA_iam['poa_back_sky_diffuse_horizon']
    
    # IAM for ground diffuse radiation components
    irradPOA_iam['poa_front_ground_diffuse'] = irradPOA_near_shade['poa_front_ground_diffuse'] * iam_diffuse_front_ground
    irradPOA_iam['poa_back_ground_diffuse'] = irradPOA_near_shade['poa_back_ground_diffuse'] * iam_diffuse_back_ground
    
    
    # complete the irradPOA_near_shade variable with remaining components.
    irradPOA_iam['poa_front_diffuse'] = irradPOA_iam['poa_front_sky_diffuse'] + irradPOA_iam['poa_front_ground_diffuse']
    irradPOA_iam['poa_back_diffuse'] = irradPOA_iam['poa_back_sky_diffuse'] + irradPOA_iam['poa_back_ground_diffuse']
    irradPOA_iam['poa_front'] = irradPOA_iam['poa_front_direct'] + irradPOA_iam['poa_front_diffuse']
    irradPOA_iam['poa_back'] = irradPOA_iam['poa_back_direct'] + irradPOA_iam['poa_back_diffuse']
    irradPOA_iam['poa_global'] = irradPOA_iam['poa_front'] + irradPOA_iam['poa_back']
    
    
    # To respect the PVsyst convention, the IAM loss is computed respect to the front side only.
    # However, irradPOA_iam already contains the components for the back side including the shading loss and IAM. Those components will be used to calculate the POA gain on the back side later
    RESULTS['IAM loss'] = np.nansum(irradPOA_iam['poa_front'])/np.nansum(irradPOA_near_shade['poa_front'])-1
    
    
    
    # =================================================================================================
    # 16 - SOILING LOSS AND POA GAIN BACK SIDE
    # =================================================================================================
    
    # effective irradiance is the irradiance hitting the PV module after IAM, soiling and snow loss
    # is the radiation that can be converted to energy
    # snow loss show be added here in the future
    
    # multiply each item in irradPOA_iam by the soiling losses
    irradPOA_soiling = (1+RESULTS['Soiling loss']) * irradPOA_iam
    
    # NOTE: the soiling loss is a uniform factor, so its percentage does not change. It is however applied to
    # the front alone in the loss chain, for the same reason explained at 'IAM loss' above
    
    # the gain of back side are finally calculated
    # NOTE: this is the OPTICAL gain of the rear side, before the bifaciality factor, so that it can be
    # compared one to one with the 'Global Irradiance on rear side' line of a PVsyst loss diagram.
    # The bifaciality factor is an optical to electrical conversion coefficient, applied further below together
    # with the STC efficiency, exactly as PVsyst does: it is therefore NOT part of the loss chain
    RESULTS['POA gain back'] = np.nansum(irradPOA_soiling['poa_back'])/np.nansum(irradPOA_soiling['poa_front'])
    
    
    
    """ SNOW LOSSES - TO BE IMPLEMENTED IN THE FUTURE
    if latitude <=55:
        snow_loss = 0     # soiling loss expresses in % - default is 2
    else:
        snow_loss = max(3, (latitude-55)/5*3)
    """
    
    
    # =================================================================================================
    # 17 - EFFECTIVE IRRADIANCE
    # =================================================================================================
    
    # apply bifaciality to back irradiance. All front components of irradEffective are the same as irradPOA_soiling. Bifaciality is the ratio between STF efficiency of front side and back side
    irradPOA_effective = pd.DataFrame(columns = POA_columns)
    irradPOA_effective['poa_front_direct'] = irradPOA_soiling['poa_front_direct']
    irradPOA_effective['poa_front_sky_diffuse_isotropic'] = irradPOA_soiling['poa_front_sky_diffuse_isotropic'] 
    irradPOA_effective['poa_front_sky_diffuse_circumsolar'] = irradPOA_soiling['poa_front_sky_diffuse_circumsolar'] 
    irradPOA_effective['poa_front_sky_diffuse_horizon'] = irradPOA_soiling['poa_front_sky_diffuse_horizon']
    irradPOA_effective['poa_front_sky_diffuse'] = irradPOA_soiling['poa_front_sky_diffuse']
    irradPOA_effective['poa_front_ground_diffuse'] = irradPOA_soiling['poa_front_ground_diffuse']
    irradPOA_effective['poa_front_diffuse'] = irradPOA_soiling['poa_front_diffuse']
    irradPOA_effective['poa_front'] = irradPOA_soiling['poa_front']
    
    irradPOA_effective['poa_back_direct'] = irradPOA_soiling['poa_back_direct'] * bifaciality 
    irradPOA_effective['poa_back_sky_diffuse_isotropic'] = irradPOA_soiling['poa_back_sky_diffuse_isotropic'] * bifaciality 
    irradPOA_effective['poa_back_sky_diffuse_circumsolar'] = irradPOA_soiling['poa_back_sky_diffuse_circumsolar'] * bifaciality 
    irradPOA_effective['poa_back_sky_diffuse_horizon'] = irradPOA_soiling['poa_back_sky_diffuse_horizon'] * bifaciality
    irradPOA_effective['poa_back_sky_diffuse'] = irradPOA_soiling['poa_back_sky_diffuse'] * bifaciality
    irradPOA_effective['poa_back_ground_diffuse'] = irradPOA_soiling['poa_back_ground_diffuse'] * bifaciality
    irradPOA_effective['poa_back_diffuse'] = irradPOA_soiling['poa_back_diffuse'] * bifaciality
    irradPOA_effective['poa_back'] =  irradPOA_soiling['poa_back'] * bifaciality
    
    irradPOA_effective['poa_global'] = irradPOA_effective['poa_front'] + irradPOA_effective['poa_back']
   
    
   
    # ================================================================================================= 
    # 18 - IRRADIANCE AND TEMPERATURE LOSS - PVSYST METHOD
    # =================================================================================================

    # this section computes the loss due to operating at non-STC condition, for both temperature and irradiation, as done on PVsyst
    
    # 18.1 - REAL MODULE CONDITIONS
    #first we compute the temperature of the cell using PVsyst method
    temp_cell_columns = ['effetive',
                         'unshaded',
                         'shaded',
                         ]
    
    temp_cell = pd.DataFrame(columns = temp_cell_columns)
    
    temp_cell['effetive'] = temperature.pvsyst_cell(poa_global = irradPOA_effective['poa_global'],
                                                    temp_air = temp_air,
                                                    wind_speed = wind_speed, 
                                                    u_c = thermal_loss_factor,            # Combined heat loss factor coefficient. The default value of 29 is representative of freestanding modules with the rear surfaces exposed to open air (e.g., rack mounted) 
                                                    u_v = 0.0,                            # Combined heat loss factor influenced by wind
                                                    module_efficiency = 0.1,              # module external efficiency. Is computed as Pdc/(poa_global_irradiance * module_area), however Pdc is not know a priori... PVsyst default value is 0.1
                                                    alpha_absorption = alpha_absorption   # PV module absorption coefficient = 1-reflectivity
                                                    )
    
    # then, it is calculated the parameters for single diode equations
    pvsyst_params = pvsystem.calcparams_pvsyst(effective_irradiance = irradPOA_effective['poa_global'],
                                               temp_cell = temp_cell['effetive'],  # temperature of PV cell using PV syst method
                                               alpha_sc = alpha_sc,                # The short-circuit current temperature coefficient of the module in units of A/C (ampere/Celsius). On PAN is reported as 'muIsc'. For Trina NEG 685 is 6.57 A/C
                                               gamma_ref = gamma_ref,              # The diode ideality factor. 1.02 is stadard pvsyst model            
                                               mu_gamma = mu_gamma,                # The temperature coefficient for the diode ideality factor, 1/K. For Trina NEG 685 is -0.0002 1/K
                                               I_L_ref = I_L_ref,                  # The light-generated current (or photocurrent) at reference conditions, equivalent to short circuit current, in amperes. for Trina NEG 685 is 18.25 A 
                                               I_o_ref = I_o_ref,                  # The dark or diode reverse saturation current at reference conditions, in amperes. 0.018 nA is standard PVsyst value
                                               R_sh_ref = R_sh_ref,                # The shunt resistance at reference conditions, in ohms. 800 is stadard PVsyst assumption
                                               R_sh_0 = R_sh_0,                    # The shunt resistance at zero irradiance conditions, in ohms. 80000 is stadard PVsyst assumption
                                               R_s = R_s,                          # The series resistance at reference conditions, in ohms. 0.168 is stadard value for Trina NEG 685
                                               cells_in_series = cells_in_series,  # The number of cells connected in series. This function assume full-cut cells. in halfcut module, cells are half the size and so actual number of cell is double respect to full-cut. In this function, alway putt the full-cut equivalent number of cell (in case of half-cut, take actual number and divide by 2)
                                               R_sh_exp = R_sh_exp,                # The exponent in the equation for shunt resistance, unitless. Defaults to 5.5
                                               EgRef = EgRef,                      # The energy bandgap at reference temperature in units of eV. 1.121 eV for crystalline silicon. EgRef must be >0
                                               irrad_ref = irrad_ref,              # STC irradiance in W/m^2
                                               temp_ref = temp_ref,                # STC temperature in °C
                                               )

    
    # function pvsystem.calcparams_pvsyst above return RESULTS in a tuple. RESULTS are rearranged in a Dataframe for convenience use
    # NaN values are replaced with 0 to allow solver for I-V curve to work
    iv_params = pd.DataFrame()
    iv_params['photocurrent'] = pvsyst_params[0]
    iv_params['saturation_current'] = pvsyst_params[1]
    iv_params['resistance_series'] = pvsyst_params[2]
    iv_params['resistance_shunt'] = pvsyst_params[3]
    iv_params['nNsVth'] = pvsyst_params[4]
    
    # the actual maximum power point of the module at the module temperature is calculated
    max_power_point = pvsystem.max_power_point(photocurrent = iv_params['photocurrent'],
                                               saturation_current = iv_params['saturation_current'], 
                                               resistance_series = iv_params['resistance_series'], 
                                               resistance_shunt = iv_params['resistance_shunt'], 
                                               nNsVth = iv_params['nNsVth'],
                                               method ='chandrupatla')     # either 'chandrupatla','brentq' or 'newton'. 'chandrupatla' slow and guarantee to converge, brent’s method mid-speed and guarantee to converge (not working here...). 'newton' is fast but NOT guarantee to converge
    
    module_real_power = max_power_point['p_mp']
    
    #18.2 - NON IDEALITY LOSS CALCULATION
    # The Theoretical Power power is calculated (STC Efficiency * G_effective)
    theoretical_power = irradPOA_effective['poa_global'] * mod_height * mod_width * stc_eff

    # RESULTS['Irradiance loss'] = np.nansum(module_ideal_power_25) / np.nansum(theoretical_power) - 1
    # RESULTS['Temperature loss'] = np.nansum(module_real_power) / np.nansum(module_ideal_power_25) - 1
    RESULTS['Temp + Irr loss'] = np.nansum(module_real_power) / np.nansum(theoretical_power) - 1   
    
    
    # ================================================================================================= 
    # 19 - PARTIAL SHADING - ELECTRICAL EFFECTS LOSS - EMPIRIC CALCULATION
    # =================================================================================================
    # Electrical effects happens when there is a partial shading of a substring, causing irradiation mismatch and hotspots.
    # Substring that are fully shaded have been accounted already in the shading loss.
    # Here we account for those substring that are partially shaded and partially unshaded
    
    # In this method, the partial shading loss is caluclated proportional to the area of the module shaded
    
    #calculate number of substrings and the width in  [m] of a module cell (in the direction of pvrow_width)
    if is_portrait and is_halfcut:        # P oriented halfcut module
        total_substrings = cells_in_parall*n_mod
        cell_width = mod_height/(cells_in_series/substrings_per_module)
        
    elif is_portrait and not is_halfcut:  # P oriented fullcut module
        total_substrings = n_mod
        cell_width = mod_height/(cells_in_series/substrings_per_module)
        
    elif not is_portrait:  # L oriented module, both halfcut and fullcut
        total_substrings = substrings_per_module*n_mod
        cell_width = mod_width/(2*substrings_per_module)
            
    # calculate the shaded fraction of the partially shaded substring (take only the decimal part of the total shaded area)
    extra_shaded_substring = SF_direct_front*total_substrings-np.floor(SF_direct_front*total_substrings)
    
    #calculate the size of one cell respect to the sixe of one substring
    cell_relative_width = cell_width/(pvrow_width/total_substrings)
    
    # calculate power loss due to electrical effect:
    # if the shaded part is less than the lenght a cell, then the substring power loss is proportional to the shaded area of the cell
    # if the shaded part is more than the lenght a cell, then all direct component is lost
    
    power_loss_partial_shading = np.where(extra_shaded_substring<cell_relative_width,
                                          extra_shaded_substring/cell_relative_width/total_substrings*(irradPOA_effective['poa_front_direct'] + irradPOA_effective['poa_front_sky_diffuse_circumsolar'])/irradPOA_effective['poa_global'],
                                          1/total_substrings*(irradPOA_effective['poa_front_direct'] + irradPOA_effective['poa_front_sky_diffuse_circumsolar'])/irradPOA_effective['poa_global']
                                          )
    
    
    RESULTS['Partial shading loss'] = np.nansum(power_loss_partial_shading * module_real_power)/np.nansum(module_real_power)*(-1)
    
    # compute the total loss inside PV module (electrical effect, module quality, LID, mismatch front, mismatch back)
    module_internal_loss = (1-power_loss_partial_shading) * (1-module_quality/100) * (1-lid_loss/100) * (1-mismatch_loss_front/100) * (1-mismatch_loss_back/100)
    
    
    
    # =================================================================================================
    # 20 - STRINGING
    # =================================================================================================    
    """
    PVyield uses a basic stringing mechanism that consist of 2xinverter configuration, one 'big' inverter with more PV strings and one 'small' inverter with less PV string.
    An inverter configuration is defined by:
        - the number of inverter of that configuration
        - the number of string connected to that inverter
        - the MPPT configuration for that inverter
    All inverter configuration infos above are included in the class InverterConfig
    
    The number of inverter with more PV strings and less PV strings are calculated to match the total number of string in the project.
    E.g. if the project has 107 PV strings to split among 12 inverter, the simulation will assume:
        'big' inverter configuration:   11x inverter of 9x string each
        'small' inverter configuration:  1x inverter of 8x string each
        Infact: 11x9 + 1x8 = 107
    NOTE: the difference in strings per inverter between the 'big' and the 'small' configuration is always 1
    
    Regarding the MPPT configuration of the inverters, it is assumed each inverter MPPT input can accept either 'less' or 'more' strings as input.
    'less' is the number of string that if applied to each MPPT that will not be enough to achieve the DC/AC ratio
    'more' is the number of string that if applied to each MPPT that will exceed the DC/AC ratio
    
    E.g. if an inverter has 12 MPPT inputs each with 3 string max and I need 32 string per inverter to achieve my DC/AC ratio:
    Therefore, I have a limitation of 36 string input on the inverter, and I need to connect 32 string
        'less' =   4 MPPT inputs of 2 string
        'more' =   8 MPPT inputs of 3 strings
    infact: 4*2 + 8*3 = 32
        
    the code will RESULTS in:
    n_mppt_less_string = 4   n_less_string = 2
    n_mppt_more_string = 8   n_more_string = 3
    
    NOTE: n_more_input is always n_less_input + 1
    """
    
    # caluclate number of inverters in project and total string in project. valuses are rounded to next integer to allow for divmod to work
    n_inv_plant = int(np.ceil(p_ac_connection/p_ac_inv_nom))    # number of plant inverters to obtain enough grid connection power. 
    tot_string_plant = int(round(p_dc_peak/modules_per_string/Pnom))  # total numer of PV strings in project, rounded to next integer
    
    if tot_string_plant < n_inv_plant:
        PVsim_err_msg += 'ERROR! The plant holds less PV strings than inverters, so at least one inverter would carry no string at all. Increase the DC power or reduce the inverter power\n\n'

    # n_string_base is the number of strings every inverter carries, n_extra_string is how many carry one more
    n_string_base, n_extra_string = divmod(tot_string_plant, n_inv_plant)
    
    InverterConfig = {'big': {},
                      'small': {},
                      }
    
    if n_extra_string == 0:
        # the strings divide evenly: every inverter carries the same number of strings and there is no 'small'
        # configuration at all. Its number of inverters is set to zero, so that it carries no weight
        InverterConfig['big']['n_string'] = n_string_base
        InverterConfig['big']['n_inv'] = n_inv_plant
        InverterConfig['small']['n_string'] = n_string_base
        InverterConfig['small']['n_inv'] = 0
    else:
        # n_extra_string inverters carry one string more than the others
        InverterConfig['big']['n_string'] = n_string_base + 1
        InverterConfig['big']['n_inv'] = n_extra_string
        InverterConfig['small']['n_string'] = n_string_base
        InverterConfig['small']['n_inv'] = n_inv_plant - n_extra_string
    
    # safety net: the strings distributed among the inverters must add up to the strings of the plant
    distributed_strings = (InverterConfig['big']['n_string'] * InverterConfig['big']['n_inv']
                           + InverterConfig['small']['n_string'] * InverterConfig['small']['n_inv'])
    if distributed_strings != tot_string_plant:
        PVsim_err_msg += 'ERROR! Internal error while distributing the PV strings among the inverters\n\n'
    
    # distribute the strings of each inverter among its MPPT inputs, following the same logic
    for config in ['big', 'small']:
        n_less_string, n_mppt_more_string = divmod(int(InverterConfig[config]['n_string']), int(n_mppt))
        
        InverterConfig[config]['n_less_string'] = n_less_string                      # number of string connected to MPPT with less string connected
        InverterConfig[config]['n_more_string'] = n_less_string + 1                  # number of string connected to MPPT with more string connected
        InverterConfig[config]['n_mppt_more_string'] = n_mppt_more_string            # number of MPPT input with the max amout of string connected
        InverterConfig[config]['n_mppt_less_string'] = n_mppt - n_mppt_more_string   # number of MPPT input with the minimum amout of string connected
    
    
    
    # =================================================================================================
    # 21 - INVERTER MPPT FINDING
    # =================================================================================================
    
    # 21.1 - PREALLOCATE VARIABLES, 
    inverter_columns = ['vmp_string',               # [V] the string voltage at module max power point conditions (=module max power point voltage * module per string), before applying MPPT limitations
                        'imp_string',               # [A] the module current at module max power point conditions, before applying MPPT limitations
                        'pmp_module',               # [W] module power at module max power point conditions, before applying MPPT limitations
                        'i_max_mppt',                # [A] max MPPT current limit applicable to the string (MPPT current limit divided by number of strings in that MPPT)
                        'vmp_mppt_i_max_limit',     # [V] inverter MPPT input lower voltage threeshold due to limitation in maximum current per each MPPT inverter input
                        'vmp_mppt_i_inv_limit',     # [V] inverter MPPT input lower threshold limit as per inverter specifications
                        'vmp_mppt_low_limit',       # [V] actual inverter MPPT input lower threshold limit , calculated as the maximum between vmp_mppt_i_max_limit and vmp_mppt_low_limit
                        'vmp_mppt_high_limit',      # [V] inverter MPPT input upeer threshold limit as per inverter specifications 
                        'vmp_mppt',                 # [V] inverter MPPT input voltage, after all voltage contraints have been accounted for
                        'imp_mppt_ideal',           # [A] MPPT DC current each string, before accounting internal PV losses (LID, quality loss, mismatch front, mismatch back, partial shading)
                        'imp_mppt',                 # [A] MPPT DC current each string, after accounting internal PV losses (LID, quality loss, mismatch front, mismatch back, partial shading)
                        'pmp_mppt_module',          # [W] power of PV module at MPPT inverter condition
                        'dc_ohm_powr_loss',         # [Wh] DC ohmic loss in the strings in Wh
                        'p_dc_string_out',          # [W] DC power at the terminals of the PV modules, before the string cable loss
                        'dc_ohm_loss',              # [-] ohmic loss in the string in [0-1] range
                        'p_dc_inv_in',              # [W] inverter DC input power
                        'p_ac_inv_out',             # [W] inverter AC output power
                        'p_dc_inv_in_ideal',        # [W] inverter DC input power for an unconstrianed inverter
                        'p_ac_inv_out',             # [W] inverter AC output power for an unconstrianed inverter
                        ]        
    
    InverterConfig['big']['mppt_less_input'] = pd.DataFrame(index = local_time, columns = inverter_columns)
    InverterConfig['big']['mppt_more_input'] = pd.DataFrame(index = local_time, columns = inverter_columns)
    InverterConfig['small']['mppt_less_input'] = pd.DataFrame(index = local_time, columns = inverter_columns)
    InverterConfig['small']['mppt_more_input'] = pd.DataFrame(index = local_time, columns = inverter_columns)
    
    
    # preallocate the MPPT voltage of the string
    InverterConfig['big']['mppt_less_input']['vmp_string'] = max_power_point ['v_mp'] * modules_per_string
    InverterConfig['big']['mppt_more_input']['vmp_string'] = max_power_point ['v_mp'] * modules_per_string
    InverterConfig['big']['mppt_less_input']['imp_string'] = max_power_point ['i_mp'] * module_internal_loss
    InverterConfig['big']['mppt_more_input']['imp_string'] = max_power_point ['i_mp'] * module_internal_loss
    InverterConfig['big']['mppt_less_input']['pmp_module'] = max_power_point ['p_mp'] * module_internal_loss
    InverterConfig['big']['mppt_more_input']['pmp_module'] = max_power_point ['p_mp'] * module_internal_loss  
    InverterConfig['big']['mppt_less_input']['vmp_mppt_i_inv_limit'] = v_min_mppt
    InverterConfig['big']['mppt_more_input']['vmp_mppt_i_inv_limit'] = v_min_mppt
    InverterConfig['big']['mppt_less_input']['vmp_mppt_high_limit'] = v_max_mppt
    InverterConfig['big']['mppt_more_input']['vmp_mppt_high_limit'] = v_max_mppt
    InverterConfig['big']['p_dc_inv_in_ideal'] = InverterConfig['big']['mppt_more_input']['pmp_module'] * modules_per_string * InverterConfig['big']['n_string']
    

    InverterConfig['small']['mppt_less_input']['vmp_string'] = max_power_point ['v_mp'] * modules_per_string
    InverterConfig['small']['mppt_more_input']['vmp_string'] = max_power_point ['v_mp'] * modules_per_string
    InverterConfig['small']['mppt_less_input']['imp_string'] = max_power_point ['i_mp'] * module_internal_loss
    InverterConfig['small']['mppt_more_input']['imp_string'] = max_power_point ['i_mp'] * module_internal_loss
    InverterConfig['small']['mppt_less_input']['pmp_module'] = max_power_point ['p_mp'] * module_internal_loss
    InverterConfig['small']['mppt_more_input']['pmp_module'] = max_power_point ['p_mp'] * module_internal_loss 
    InverterConfig['small']['mppt_less_input']['vmp_mppt_i_inv_limit'] = v_min_mppt
    InverterConfig['small']['mppt_more_input']['vmp_mppt_i_inv_limit'] = v_min_mppt
    InverterConfig['small']['mppt_less_input']['vmp_mppt_high_limit'] = v_max_mppt
    InverterConfig['small']['mppt_more_input']['vmp_mppt_high_limit'] = v_max_mppt
    InverterConfig['small']['p_dc_inv_in_ideal'] = InverterConfig['small']['mppt_more_input']['pmp_module'] * modules_per_string * InverterConfig['small']['n_string']    
    
    # 21.2 - COMPUTE MPPT VOLTAGE BOUNDARIES
    """
    The inverter MPPT input has a maximum curent limit. This limit is for the whole MPPT.
    So if there are multiple strings connected to that MPPT input, it each string can carry a maximum value equal to the inverter MPPT current limit divided by the string connected to the MPPT.
    e.g. if the inverter MPPT input has a limit of 30A, and can connect 3 strings max, each string can carry 10A max
    
    We need to be sure that the MPP of the module generates a current that is lower than this limit (10A in example abvove).
    Therefore, the inverter MPPT voltage scan cannot select a voltage corresponding to a current higher than this limit.
    Since the IV curve is monotonically decreasing, this max current limit generate a lower voltage limit for the inverter MPPT scan.
    
    In other words, the inverter MPPT max curent limit generate a MPPT voltage lower limit. If the voltage of the string is lower than this limit, the current produced will be higher than the acceptable inverter MPPT current limit
    The voltage limit due to current is compared with the voltage limit of the MPPT, and the most stringent of the two (so the maximum among the two) is taken in to account as MPPT voltage scan limit
    """
    
    if not(InverterConfig['big']['n_less_string']==0):
        # starting from the current limit of the inverter MPPT, compute the max MPPT current limit applicable to the string, for the MMPT input with less string
        InverterConfig['big']['mppt_less_input']['i_max_mppt'] = i_max_mppt/module_internal_loss/InverterConfig['big']['n_less_string']
        
        # compute lower boundary of MPPT voltage due to max current limit, for MMPT input with less string
        # if the max current limit exceed the short circuit current of the module, then the max current limit is set to the short circuit current of the module
        # by setting the max current limit to the short circuit, the v_from_i function will return a voltage equal to 0
        InverterConfig['big']['mppt_less_input']['i_max_mppt'] = np.where(InverterConfig['big']['mppt_less_input']['i_max_mppt'] >= iv_params['photocurrent'], iv_params['photocurrent'], InverterConfig['big']['mppt_less_input']['i_max_mppt'])
        InverterConfig['big']['mppt_less_input']['vmp_mppt_i_max_limit']  = modules_per_string *pvsystem.v_from_i(current = InverterConfig['big']['mppt_less_input']['i_max_mppt'],
                                                                                                                  photocurrent = iv_params['photocurrent'],
                                                                                                                  saturation_current = iv_params['saturation_current'], 
                                                                                                                  resistance_series = iv_params['resistance_series'], 
                                                                                                                  resistance_shunt = iv_params['resistance_shunt'], 
                                                                                                                  nNsVth = iv_params['nNsVth'],
                                                                                                                  method='lambertw',      # Method to use: 'lambertw', 'chandrupatla', or 'brentq'. Note: 'brentq' is limited to 1st quadrant only
                                                                                                                  )
    
    else:
        InverterConfig['big']['mppt_less_input']['vmp_mppt_i_max_limit'] = 0
    
    
    if not(InverterConfig['small']['n_less_string']==0):
        # starting from the current limit of the inverter MPPT, compute the max MPPT current limit applicable to the string, for the MMPT input with less string
        InverterConfig['small']['mppt_less_input']['i_max_mppt'] = i_max_mppt/module_internal_loss/InverterConfig['small']['n_less_string']
        
        # compute lower boundary of MPPT voltage due to max current limit, for MMPT input with less string
        # if the max current limit exceed the short circuit current of the module, then the max current limit is set to the short circuit current of the module
        # by setting the max current limit to the short circuit, the v_from_i function will return a voltage equal to 0
        InverterConfig['small']['mppt_less_input']['i_max_mppt'] = np.where(InverterConfig['small']['mppt_less_input']['i_max_mppt'] >= iv_params['photocurrent'], iv_params['photocurrent'], InverterConfig['small']['mppt_less_input']['i_max_mppt'])
        InverterConfig['small']['mppt_less_input']['vmp_mppt_i_max_limit']  = modules_per_string *pvsystem.v_from_i(current = InverterConfig['small']['mppt_less_input']['i_max_mppt'],
                                                                                                                    photocurrent = iv_params['photocurrent'],
                                                                                                                    saturation_current = iv_params['saturation_current'], 
                                                                                                                    resistance_series = iv_params['resistance_series'], 
                                                                                                                    resistance_shunt = iv_params['resistance_shunt'], 
                                                                                                                    nNsVth = iv_params['nNsVth'],
                                                                                                                    method='lambertw',      # Method to use: 'lambertw', 'newton', or 'brentq'. Note: 'brentq' is limited to 1st quadrant only
                                                                                                                    )
    
    else:
        InverterConfig['small']['mppt_less_input']['vmp_mppt_i_max_limit'] = 0
    
    
    # starting from the current limit of the inverter MPPT, compute the max MPPT current limit applicable to the string, for the MMPT input with more string
    InverterConfig['big']['mppt_more_input']['i_max_mppt'] = i_max_mppt/module_internal_loss/InverterConfig['big']['n_more_string']
    InverterConfig['small']['mppt_more_input']['i_max_mppt'] = i_max_mppt/module_internal_loss/InverterConfig['small']['n_more_string']
    
    # compute lower boundary of MPPT voltage due to max current limit, for MMPT input with more string
    # if the max current limit exceed the short circuit current of the module, then the max current limit is set to the short circuit current of the module
    # by setting the max current limit to the short circuit, the v_from_i function will return a voltage equal to 0
    InverterConfig['big']['mppt_more_input']['i_max_mppt'] = np.where(InverterConfig['big']['mppt_more_input']['i_max_mppt'] >= iv_params['photocurrent'], iv_params['photocurrent'], InverterConfig['big']['mppt_more_input']['i_max_mppt'])
    InverterConfig['big']['mppt_more_input']['vmp_mppt_i_max_limit'] = modules_per_string * pvsystem.v_from_i(current = InverterConfig['big']['mppt_more_input']['i_max_mppt'],
                                                                                                              photocurrent = iv_params['photocurrent'],
                                                                                                              saturation_current = iv_params['saturation_current'], 
                                                                                                              resistance_series = iv_params['resistance_series'], 
                                                                                                              resistance_shunt = iv_params['resistance_shunt'], 
                                                                                                              nNsVth = iv_params['nNsVth'],
                                                                                                              method='lambertw',      # Method to use: 'lambertw', 'newton', or 'brentq'. Note: 'brentq' is limited to 1st quadrant only
                                                                                                              )
    
    InverterConfig['small']['mppt_more_input']['i_max_mppt'] = np.where(InverterConfig['small']['mppt_more_input']['i_max_mppt'] >= iv_params['photocurrent'], iv_params['photocurrent'], InverterConfig['small']['mppt_more_input']['i_max_mppt'])
    InverterConfig['small']['mppt_more_input']['vmp_mppt_i_max_limit'] = modules_per_string * pvsystem.v_from_i(current = InverterConfig['small']['mppt_more_input']['i_max_mppt'],
                                                                                                                photocurrent = iv_params['photocurrent'],
                                                                                                                saturation_current = iv_params['saturation_current'], 
                                                                                                                resistance_series = iv_params['resistance_series'], 
                                                                                                                resistance_shunt = iv_params['resistance_shunt'], 
                                                                                                                nNsVth = iv_params['nNsVth'],
                                                                                                                method='lambertw',      # Method to use: 'lambertw', 'newton', or 'brentq'. Note: 'brentq' is limited to 1st quadrant only
                                                                                                                )
    
    # 21.3 - RECOMPUTE MPPT VOLTAGE TO MAKE SURE IT IS INSIDE BOUNDARIES
    # compute the new minimum voltage MPPT limit by taking the maximum between the limit derived from the max current limit and the one from inverter MPPT specification
    InverterConfig['big']['mppt_less_input']['vmp_mppt_low_limit'] = np.where(InverterConfig['big']['mppt_less_input']['vmp_mppt_i_max_limit'] > v_min_mppt, InverterConfig['big']['mppt_less_input']['vmp_mppt_i_max_limit'], v_min_mppt)
    InverterConfig['big']['mppt_more_input']['vmp_mppt_low_limit'] = np.where(InverterConfig['big']['mppt_more_input']['vmp_mppt_i_max_limit'] > v_min_mppt, InverterConfig['big']['mppt_more_input']['vmp_mppt_i_max_limit'], v_min_mppt)
    
    InverterConfig['small']['mppt_less_input']['vmp_mppt_low_limit'] = np.where(InverterConfig['small']['mppt_less_input']['vmp_mppt_i_max_limit'] > v_min_mppt, InverterConfig['small']['mppt_less_input']['vmp_mppt_i_max_limit'], v_min_mppt)
    InverterConfig['small']['mppt_more_input']['vmp_mppt_low_limit'] = np.where(InverterConfig['small']['mppt_more_input']['vmp_mppt_i_max_limit'] > v_min_mppt, InverterConfig['small']['mppt_more_input']['vmp_mppt_i_max_limit'], v_min_mppt)
    
    # compute the actualg MPP voltage by clipping the original string value to the lower and upper MPPT voltage boundaries
    InverterConfig['big']['mppt_less_input']['vmp_mppt'] = np.clip(InverterConfig['big']['mppt_less_input']['vmp_string'], InverterConfig['big']['mppt_less_input']['vmp_mppt_low_limit'], InverterConfig['big']['mppt_less_input']['vmp_mppt_high_limit'])  
    InverterConfig['big']['mppt_more_input']['vmp_mppt'] = np.clip(InverterConfig['big']['mppt_more_input']['vmp_string'], InverterConfig['big']['mppt_more_input']['vmp_mppt_low_limit'], InverterConfig['big']['mppt_more_input']['vmp_mppt_high_limit']) 
    
    InverterConfig['small']['mppt_less_input']['vmp_mppt'] = np.clip(InverterConfig['small']['mppt_less_input']['vmp_string'], InverterConfig['small']['mppt_less_input']['vmp_mppt_low_limit'], InverterConfig['small']['mppt_less_input']['vmp_mppt_high_limit'])  
    InverterConfig['small']['mppt_more_input']['vmp_mppt'] = np.clip(InverterConfig['small']['mppt_more_input']['vmp_string'], InverterConfig['small']['mppt_more_input']['vmp_mppt_low_limit'], InverterConfig['small']['mppt_more_input']['vmp_mppt_high_limit']) 
    
    
    # 21.4 - COMPUTE MPPT CURRENT FROM VOLTAGE
    # calculate the current generated by the module at the MPP voltage found
    InverterConfig['big']['mppt_less_input']['imp_mppt_ideal'] = pvsystem.i_from_v(voltage = InverterConfig['big']['mppt_less_input']['vmp_mppt']/modules_per_string,
                                                                                   photocurrent = iv_params['photocurrent'],
                                                                                   saturation_current = iv_params['saturation_current'], 
                                                                                   resistance_series = iv_params['resistance_series'], 
                                                                                   resistance_shunt = iv_params['resistance_shunt'], 
                                                                                   nNsVth = iv_params['nNsVth'],
                                                                                   method='lambertw',      # Method to use: 'lambertw', 'newton', or 'brentq'. Note: 'brentq' is limited to 1st quadrant only
                                                                                   )
    
    InverterConfig['big']['mppt_more_input']['imp_mppt_ideal'] = pvsystem.i_from_v(voltage = InverterConfig['big']['mppt_more_input']['vmp_mppt']/modules_per_string,
                                                                                   photocurrent = iv_params['photocurrent'],
                                                                                   saturation_current = iv_params['saturation_current'], 
                                                                                   resistance_series = iv_params['resistance_series'], 
                                                                                   resistance_shunt = iv_params['resistance_shunt'], 
                                                                                   nNsVth = iv_params['nNsVth'],
                                                                                   method='lambertw',      # Method to use: 'lambertw', 'newton', or 'brentq'. Note: 'brentq' is limited to 1st quadrant only
                                                                                   )
    
    InverterConfig['small']['mppt_less_input']['imp_mppt_ideal'] = pvsystem.i_from_v(voltage = InverterConfig['small']['mppt_less_input']['vmp_mppt']/modules_per_string,
                                                                                     photocurrent = iv_params['photocurrent'],
                                                                                     saturation_current = iv_params['saturation_current'], 
                                                                                     resistance_series = iv_params['resistance_series'], 
                                                                                     resistance_shunt = iv_params['resistance_shunt'], 
                                                                                     nNsVth = iv_params['nNsVth'],
                                                                                     method='lambertw',      # Method to use: 'lambertw', 'newton', or 'brentq'. Note: 'brentq' is limited to 1st quadrant only
                                                                                     )
    
    InverterConfig['small']['mppt_more_input']['imp_mppt_ideal'] = pvsystem.i_from_v(voltage = InverterConfig['small']['mppt_more_input']['vmp_mppt']/modules_per_string,
                                                                                     photocurrent = iv_params['photocurrent'],
                                                                                     saturation_current = iv_params['saturation_current'], 
                                                                                     resistance_series = iv_params['resistance_series'], 
                                                                                     resistance_shunt = iv_params['resistance_shunt'], 
                                                                                     nNsVth = iv_params['nNsVth'],
                                                                                     method='lambertw',      # Method to use: 'lambertw', 'newton', or 'brentq'. Note: 'brentq' is limited to 1st quadrant only
                                                                                     )
    
    
    # the actual current flowing in the string is lower than the resul from IV curve, due to module internal loss (lID, quality loss, front mismatch, back mismatch, partial shading)
    # the MPP current found above is hence reduced by a factor equal to the PM module internal losses
    InverterConfig['big']['mppt_less_input']['imp_mppt'] = InverterConfig['big']['mppt_less_input']['imp_mppt_ideal'] * module_internal_loss
    InverterConfig['big']['mppt_more_input']['imp_mppt'] = InverterConfig['big']['mppt_more_input']['imp_mppt_ideal'] * module_internal_loss

    InverterConfig['small']['mppt_less_input']['imp_mppt'] = InverterConfig['small']['mppt_less_input']['imp_mppt_ideal'] * module_internal_loss
    InverterConfig['small']['mppt_more_input']['imp_mppt'] = InverterConfig['small']['mppt_more_input']['imp_mppt_ideal'] * module_internal_loss
    
    # the PV module power at MPPT condition ic calculated
    InverterConfig['big']['mppt_less_input']['pmp_mppt_module'] = InverterConfig['big']['mppt_less_input']['imp_mppt'] * InverterConfig['big']['mppt_less_input']['vmp_mppt'] / modules_per_string
    InverterConfig['big']['mppt_more_input']['pmp_mppt_module'] = InverterConfig['big']['mppt_more_input']['imp_mppt'] * InverterConfig['big']['mppt_more_input']['vmp_mppt'] / modules_per_string

    InverterConfig['small']['mppt_less_input']['pmp_mppt_module'] = InverterConfig['small']['mppt_less_input']['imp_mppt'] * InverterConfig['small']['mppt_less_input']['vmp_mppt'] / modules_per_string
    InverterConfig['small']['mppt_more_input']['pmp_mppt_module'] = InverterConfig['small']['mppt_more_input']['imp_mppt'] * InverterConfig['small']['mppt_more_input']['vmp_mppt'] / modules_per_string



    # =================================================================================================
    # 22 DC OHMIC LOSS
    # =================================================================================================
    
    # compute equivalent resistance needed to compute DC string cable losses
    
    InverterConfig['big']['eq_resistance_less_input'] = pvsystem.dc_ohms_from_percent(vmp_ref = vmp_ref,   # Voltage at maximum power in reference conditions [V]
                                                                                      imp_ref = imp_ref,   # Current at maximum power in reference conditions [I]
                                                                                      dc_ohmic_percent = dc_ohmic_loss_STC,
                                                                                      modules_per_string = modules_per_string, 
                                                                                      strings = InverterConfig['big']['n_less_string'],   # Number of parallel strings in the array
                                                                                      )
    
    InverterConfig['big']['eq_resistance_more_input'] = pvsystem.dc_ohms_from_percent(vmp_ref = vmp_ref,   # Voltage at maximum power in reference conditions [V]
                                                                                      imp_ref = imp_ref,   # Current at maximum power in reference conditions [I]
                                                                                      dc_ohmic_percent = dc_ohmic_loss_STC,
                                                                                      modules_per_string = modules_per_string, 
                                                                                      strings = InverterConfig['big']['n_more_string'],   # Number of parallel strings in the array
                                                                                      )
    
    
    InverterConfig['small']['eq_resistance_less_input'] = pvsystem.dc_ohms_from_percent(vmp_ref = vmp_ref,   # Voltage at maximum power in reference conditions [V]
                                                                                        imp_ref = imp_ref,   # Current at maximum power in reference conditions [I]
                                                                                        dc_ohmic_percent = dc_ohmic_loss_STC,
                                                                                        modules_per_string = modules_per_string, 
                                                                                        strings = InverterConfig['small']['n_less_string'],   # Number of parallel strings in the array
                                                                                        )
    
    InverterConfig['small']['eq_resistance_more_input'] = pvsystem.dc_ohms_from_percent(vmp_ref = vmp_ref,   # Voltage at maximum power in reference conditions [V]
                                                                                        imp_ref = imp_ref,   # Current at maximum power in reference conditions [I]
                                                                                        dc_ohmic_percent = dc_ohmic_loss_STC,
                                                                                        modules_per_string = modules_per_string, 
                                                                                        strings = InverterConfig['small']['n_more_string'],   # Number of parallel strings in the array
                                                                                        )
    
    
    # 22.1 - CALCULATE DC OHMIC LOSSES
    # NOTE: The function dc_ohmic_losses returns the DC ohmic losses in Wh
    # NOTE: the operating current is used on purpose. The STC reference is already carried by Rw, so the loss
    # correctly scales with the square of the actual current at every timestep
    for config in ['big', 'small']:
        for mppt, suffix in [('mppt_less_input', 'less'), ('mppt_more_input', 'more')]:
            n_string_input = InverterConfig[config]['n_'+suffix+'_string']
            n_mppt_input = InverterConfig[config]['n_mppt_'+suffix+'_string']
            
            if n_string_input == 0:     # no string connected to this MPPT input type, hence no DC losses
                InverterConfig[config][mppt]['dc_ohm_powr_loss'] = 0
            else:
                # current flowing through one MPPT input: the string current times the strings connected to it
                i_mppt_input = InverterConfig[config][mppt]['imp_mppt'] * n_string_input
                
                InverterConfig[config][mppt]['dc_ohm_powr_loss'] = pvsystem.dc_ohmic_losses(resistance = InverterConfig[config]['eq_resistance_'+suffix+'_input'],
                                                                                            current = i_mppt_input,
                                                                                            )
                
                # scale the loss of a single MPPT input up to all the MPPT inputs of that type
                InverterConfig[config][mppt]['dc_ohm_powr_loss'] = InverterConfig[config][mppt]['dc_ohm_powr_loss'] * n_mppt_input
            
            # power leaving the PV modules connected to this MPPT input type, at the module terminals
            InverterConfig[config][mppt]['p_dc_string_out'] = InverterConfig[config][mppt]['pmp_mppt_module'] * modules_per_string * n_string_input * n_mppt_input
            
            # NOTE ON THE OPERATING POINT: the MPPT search above does not 'feel' the cable resistance, therefore maximises the power at the MODULE terminals and ignores the cable resistance.
            # A real inverter would 'feel' the cable resistance, finding MPPT points at its own terminaly. Therefore, a real inverter MPPT would be at slightly higher voltage and drawing a slightly lower current.
            # This effect is deliberately left out: it is a second order effect accounting for a fraction (<0.01 %) of yield
            # It is possible to add this effect by adding the cable resistance to R_s before calcparams_pvsyst, which also changes the cell temperature and the MPPT search, and is not worth it for a gain of this size.

            InverterConfig[config][mppt]['p_dc_inv_in'] = (InverterConfig[config][mppt]['p_dc_string_out']
                                                           - InverterConfig[config][mppt]['dc_ohm_powr_loss'])
        
        InverterConfig[config]['dc_ohm_powr_loss'] = InverterConfig[config]['mppt_less_input']['dc_ohm_powr_loss'] + InverterConfig[config]['mppt_more_input']['dc_ohm_powr_loss']
        InverterConfig[config]['p_dc_string_out'] = InverterConfig[config]['mppt_less_input']['p_dc_string_out'] + InverterConfig[config]['mppt_more_input']['p_dc_string_out']
        InverterConfig[config]['p_dc_inv_in'] = InverterConfig[config]['mppt_less_input']['p_dc_inv_in'] + InverterConfig[config]['mppt_more_input']['p_dc_inv_in']
    
    # NOTE: dc_ohm_powr_loss_tot is no longer needed to build the loss row, the cable loss being already
    # subtracted from p_dc_inv_in above. It is kept because it is the only place where the DC cable loss
    # of the whole plant is available in Wh, which is handy when debugging
    p_dc_string_out = InverterConfig['big']['p_dc_string_out'] * InverterConfig['big']['n_inv'] + InverterConfig['small']['p_dc_string_out'] * InverterConfig['small']['n_inv']
    p_dc_inv_in = InverterConfig['big']['p_dc_inv_in'] * InverterConfig['big']['n_inv'] + InverterConfig['small']['p_dc_inv_in'] * InverterConfig['small']['n_inv'] 

    # MPPT VOLTAGE AND CURRENT CONSTRAINTS
    # The MPPT of the inverter cannot always sit on the maximum power point of the modules: its voltage window
    # and its maximum input current push the operating point away from it. The energy lost that way is the
    # difference between the power the strings deliver at the constrained operating point and the power they
    # would deliver at the true maximum power point, both taken at the module terminals.
    #
    # This loss used to be hidden inside a row called 'Inverter thresholds loss', computed as the
    # ratio between the IDEAL clipped output and the REAL pre-clipping one. Mixing the two chains gave a row
    # that did not measure its own stage, and with a real OND file it even came out POSITIVE: on the reference
    # project it reported +0.498%, i.e. a gain, where the inverter can only take energy away. Measured on that
    # same project, the ideal chain delivers 14.3397 GWh against the 14.2686 GWh of the real one, so the true
    # MPPT constraint loss is -0.496% and the clipping is nil.
    # The two effects are now separated and each row measures exactly the stage it names:
    #   'Inverter MPPT V/I loss' : computed here, on the DC side, where the constraint physically acts
    #   'Inverter clipping loss' : computed further below, on the AC side, where the output is capped
    # Both are always negative or zero, and the loss chain telescopes without needing any correction
    p_dc_string_out_ideal = InverterConfig['big']['p_dc_inv_in_ideal'] * InverterConfig['big']['n_inv'] + InverterConfig['small']['p_dc_inv_in_ideal'] * InverterConfig['small']['n_inv']
    
    # NOTE: the row is DISPLAYED further down the list, between the string cable loss and the clipping, the
    # way PVsyst orders its own loss diagram. It is COMPUTED here because this is where the two DC powers it
    # compares are available. The position in the list has no effect on the numbers: the Performance Ratio is
    # a ratio of energies, and the loss rows are referenced by name and never by position
    RESULTS['Inverter MPPT V/I loss'] = np.nansum(p_dc_string_out)/np.nansum(p_dc_string_out_ideal) - 1
    
    RESULTS['DC ohmic wiring loss'] = np.nansum(p_dc_inv_in)/np.nansum(p_dc_string_out) - 1
    
    # # compute DC ohmic losses in 0-1 range
    # InverterConfig['big']['mppt_less_input']['dc_ohm_loss'] = (-1) * InverterConfig['big']['mppt_less_input']['dc_ohm_powr_loss'] / (InverterConfig['big']['mppt_less_input']['dc_ohm_powr_loss'] + InverterConfig['big']['mppt_less_input']['p_dc_inv_in'])
    # InverterConfig['big']['mppt_more_input']['dc_ohm_loss'] = (-1) * InverterConfig['big']['mppt_more_input']['dc_ohm_powr_loss'] / (InverterConfig['big']['mppt_more_input']['dc_ohm_powr_loss'] + InverterConfig['big']['mppt_more_input']['p_dc_inv_in'])
    # InverterConfig['small']['mppt_less_input']['dc_ohm_loss'] = (-1) * InverterConfig['small']['mppt_less_input']['dc_ohm_powr_loss'] / (InverterConfig['small']['mppt_less_input']['dc_ohm_powr_loss'] + InverterConfig['small']['mppt_less_input']['p_dc_inv_in'])
    # InverterConfig['small']['mppt_more_input']['dc_ohm_loss'] = (-1) * InverterConfig['small']['mppt_more_input']['dc_ohm_powr_loss'] / (InverterConfig['small']['mppt_more_input']['dc_ohm_powr_loss'] + InverterConfig['small']['mppt_more_input']['p_dc_inv_in'])
    
    # InverterConfig['big']['dc_ohm_loss'] = (-1)*InverterConfig['big']['dc_ohm_powr_loss'] / InverterConfig['big']['p_dc_out_string_inv']
    # InverterConfig['small']['dc_ohm_loss'] = (-1)*InverterConfig['small']['dc_ohm_powr_loss'] / InverterConfig['small']['p_dc_out_string_inv']
    

    # # 22.2 CALCULATE INVERTER EFFICIENCY BASED ON Pdc_in and OUTPUT VOLTAGE
    # # compute DC power exiting the string, before DC ohmic loss
    # InverterConfig['big']['mppt_less_input']['p_dc_out_string_inv'] = InverterConfig['big']['mppt_less_input']['p_dc_inv_in'] / (1 + InverterConfig['big']['mppt_less_input']['dc_ohm_loss'])
    # InverterConfig['big']['mppt_more_input']['p_dc_out_string_inv'] = InverterConfig['big']['mppt_more_input']['p_dc_inv_in'] / (1 + InverterConfig['big']['mppt_more_input']['dc_ohm_loss'])
    # InverterConfig['small']['mppt_less_input']['p_dc_out_string_inv'] = InverterConfig['small']['mppt_less_input']['p_dc_inv_in'] / (1 + InverterConfig['small']['mppt_less_input']['dc_ohm_loss'])
    # InverterConfig['small']['mppt_more_input']['p_dc_out_string_inv'] = InverterConfig['small']['mppt_more_input']['p_dc_inv_in'] / (1 + InverterConfig['small']['mppt_more_input']['dc_ohm_loss'])
    
    # InverterConfig['big']['p_dc_out_string_inv'] = InverterConfig['big']['p_dc_inv_in'] / (1 + InverterConfig['big']['dc_ohm_loss'])
    # InverterConfig['small']['p_dc_out_string_inv'] = InverterConfig['big']['p_dc_inv_in'] /  (1 + InverterConfig['small']['dc_ohm_loss'])
    
    
    # compute output power as function of the 3xefficiency curves specified in OND: P_out = f(P_in, V_in, P_out_max, efficiency curves)
    #NOTE: the function take as input the efifciency curve expressed over the total inverter input/output power range, while 'p_dc_inv_in' is the power input for a limited set of MPPT inputs
    # therefore, p_dc_inv_in and p_ac_inv_out are scaled to assume that the whole inverter is working with the spefcific power input of the MPPT input 
    
    if INVERTER['OND_import']:
        # Efficiency calculation per MPPT input not working, it has been commented. Please leave efficiency calculation per inverter as below
        # InverterConfig['big']['mppt_less_input']['p_ac_inv_out'] = InverterConfig['big']['n_mppt_less_string'] / n_mppt * inv_Pout_from_Pin(v_eff1 = v_eff1,
        #                                                                                                                                     v_eff2 = v_eff2,
        #                                                                                                                                     v_eff3 = v_eff3,
        #                                                                                                                                     p_in_v_eff1 = p_in_v_eff1,
        #                                                                                                                                     p_in_v_eff2 = p_in_v_eff2,
        #                                                                                                                                     p_in_v_eff3 = p_in_v_eff3,
        #                                                                                                                                     p_out_v_eff1 = p_out_v_eff1,
        #                                                                                                                                     p_out_v_eff2 = p_out_v_eff2,
        #                                                                                                                                     p_out_v_eff3 = p_out_v_eff3,
        #                                                                                                                                     p_dc_inv_in = InverterConfig['big']['mppt_less_input']['p_dc_inv_in'] / InverterConfig['big']['n_mppt_less_string'] * n_mppt,
        #                                                                                                                                     v_dc_in = InverterConfig['big']['mppt_less_input']['vmp_mppt'],
        #                                                                                                                                     )
        
        # InverterConfig['big']['mppt_more_input']['p_ac_inv_out'] = InverterConfig['big']['n_mppt_more_string'] / n_mppt * inv_Pout_from_Pin(v_eff1 = v_eff1,
        #                                                                                                                                       v_eff2 = v_eff2,
        #                                                                                                                                       v_eff3 = v_eff3,
        #                                                                                                                                       p_in_v_eff1 = p_in_v_eff1,
        #                                                                                                                                       p_in_v_eff2 = p_in_v_eff2,
        #                                                                                                                                       p_in_v_eff3 = p_in_v_eff3,
        #                                                                                                                                       p_out_v_eff1 = p_out_v_eff1,
        #                                                                                                                                       p_out_v_eff2 = p_out_v_eff2,
        #                                                                                                                                       p_out_v_eff3 = p_out_v_eff3,
        #                                                                                                                                       p_dc_inv_in = InverterConfig['big']['mppt_more_input']['p_dc_inv_in'] / InverterConfig['big']['n_mppt_more_string'] * n_mppt,
        #                                                                                                                                       v_dc_in = InverterConfig['big']['mppt_more_input']['vmp_mppt'],
        #                                                                                                                                       )
        
        # InverterConfig['small']['mppt_less_input']['p_ac_inv_out'] = InverterConfig['small']['n_mppt_less_string'] / n_mppt * inv_Pout_from_Pin(v_eff1 = v_eff1,
        #                                                                                                                                         v_eff2 = v_eff2,
        #                                                                                                                                         v_eff3 = v_eff3,
        #                                                                                                                                         p_in_v_eff1 = p_in_v_eff1,
        #                                                                                                                                         p_in_v_eff2 = p_in_v_eff2,
        #                                                                                                                                         p_in_v_eff3 = p_in_v_eff3,
        #                                                                                                                                         p_out_v_eff1 = p_out_v_eff1,
        #                                                                                                                                         p_out_v_eff2 = p_out_v_eff2,
        #                                                                                                                                         p_out_v_eff3 = p_out_v_eff3,
        #                                                                                                                                         p_dc_inv_in = InverterConfig['small']['mppt_less_input']['p_dc_inv_in'] / InverterConfig['small']['n_mppt_less_string'] * n_mppt,
        #                                                                                                                                         v_dc_in = InverterConfig['small']['mppt_less_input']['vmp_mppt'],
        #                                                                                                                                         )
        
        # InverterConfig['small']['mppt_more_input']['p_ac_inv_out'] = InverterConfig['small']['n_mppt_more_string'] / n_mppt * inv_Pout_from_Pin(v_eff1 = v_eff1,
        #                                                                                                                                         v_eff2 = v_eff2,
        #                                                                                                                                         v_eff3 = v_eff3,
        #                                                                                                                                         p_in_v_eff1 = p_in_v_eff1,
        #                                                                                                                                         p_in_v_eff2 = p_in_v_eff2,
        #                                                                                                                                         p_in_v_eff3 = p_in_v_eff3,
        #                                                                                                                                         p_out_v_eff1 = p_out_v_eff1,
        #                                                                                                                                         p_out_v_eff2 = p_out_v_eff2,
        #                                                                                                                                         p_out_v_eff3 = p_out_v_eff3,
        #                                                                                                                                         p_dc_inv_in = InverterConfig['small']['mppt_more_input']['p_dc_inv_in'] / InverterConfig['small']['n_mppt_more_string'] * n_mppt,
        #                                                                                                                                         v_dc_in = InverterConfig['small']['mppt_more_input']['vmp_mppt'],
        #                                                                                                                                         )
        
        # InverterConfig['big']['p_ac_inv_out'] = InverterConfig['big']['mppt_less_input']['p_ac_inv_out'] + InverterConfig['big']['mppt_more_input']['p_ac_inv_out']
        # InverterConfig['small']['p_ac_inv_out'] = InverterConfig['small']['mppt_less_input']['p_ac_inv_out'] + InverterConfig['small']['mppt_more_input']['p_ac_inv_out']
        
        # Inverter efficiency calculation is done at inverter level, not MPPT, since MPPT calculation (commented section above) does not work well
        
        InverterConfig['big']['p_ac_inv_out'] = inv_Pout_from_Pin(v_eff1 = v_eff1,
                                                                  v_eff2 = v_eff2,
                                                                  v_eff3 = v_eff3,
                                                                  p_in_v_eff1 = p_in_v_eff1,
                                                                  p_in_v_eff2 = p_in_v_eff2,
                                                                  p_in_v_eff3 = p_in_v_eff3,
                                                                  p_out_v_eff1 = p_out_v_eff1,
                                                                  p_out_v_eff2 = p_out_v_eff2,
                                                                  p_out_v_eff3 = p_out_v_eff3,
                                                                  p_dc_inv_in = InverterConfig['big']['p_dc_inv_in'],
                                                                  v_dc_in = (InverterConfig['big']['mppt_more_input']['vmp_mppt'] * InverterConfig['big']['n_mppt_more_string'] + InverterConfig['big']['mppt_less_input']['vmp_mppt'] * InverterConfig['big']['n_mppt_less_string']) / n_mppt,
                                                                  )
        
        InverterConfig['small']['p_ac_inv_out'] = inv_Pout_from_Pin(v_eff1 = v_eff1,
                                                                    v_eff2 = v_eff2,
                                                                    v_eff3 = v_eff3,
                                                                    p_in_v_eff1 = p_in_v_eff1,
                                                                    p_in_v_eff2 = p_in_v_eff2,
                                                                    p_in_v_eff3 = p_in_v_eff3,
                                                                    p_out_v_eff1 = p_out_v_eff1,
                                                                    p_out_v_eff2 = p_out_v_eff2,
                                                                    p_out_v_eff3 = p_out_v_eff3,
                                                                    p_dc_inv_in = InverterConfig['small']['p_dc_inv_in'],
                                                                    v_dc_in = (InverterConfig['small']['mppt_more_input']['vmp_mppt'] * InverterConfig['small']['n_mppt_more_string'] + InverterConfig['small']['mppt_less_input']['vmp_mppt'] * InverterConfig['small']['n_mppt_less_string']) / n_mppt,
                                                                    )
        
        InverterConfig['big']['p_ac_inv_out_ideal'] = inv_Pout_from_Pin(v_eff1 = v_eff1,
                                                                        v_eff2 = v_eff2,
                                                                        v_eff3 = v_eff3,
                                                                        p_in_v_eff1 = p_in_v_eff1,
                                                                        p_in_v_eff2 = p_in_v_eff2,
                                                                        p_in_v_eff3 = p_in_v_eff3,
                                                                        p_out_v_eff1 = p_out_v_eff1,
                                                                        p_out_v_eff2 = p_out_v_eff2,
                                                                        p_out_v_eff3 = p_out_v_eff3,
                                                                        p_dc_inv_in = InverterConfig['big']['p_dc_inv_in_ideal'],
                                                                        v_dc_in = max_power_point ['v_mp'] * modules_per_string,
                                                                        )
        
        InverterConfig['small']['p_ac_inv_out_ideal'] = inv_Pout_from_Pin(v_eff1 = v_eff1,
                                                                    v_eff2 = v_eff2,
                                                                    v_eff3 = v_eff3,
                                                                    p_in_v_eff1 = p_in_v_eff1,
                                                                    p_in_v_eff2 = p_in_v_eff2,
                                                                    p_in_v_eff3 = p_in_v_eff3,
                                                                    p_out_v_eff1 = p_out_v_eff1,
                                                                    p_out_v_eff2 = p_out_v_eff2,
                                                                    p_out_v_eff3 = p_out_v_eff3,
                                                                    p_dc_inv_in = InverterConfig['small']['p_dc_inv_in_ideal'],
                                                                    v_dc_in = max_power_point ['v_mp'] * modules_per_string,
                                                                    )
        
        p_ac_inv_out_preclipping = InverterConfig['big']['p_ac_inv_out'] * InverterConfig['big']['n_inv'] + InverterConfig['small']['p_ac_inv_out'] * InverterConfig['small']['n_inv']
        p_ac_inv_out_preclipping_ideal = InverterConfig['big']['p_ac_inv_out_ideal'] * InverterConfig['big']['n_inv'] + InverterConfig['small']['p_ac_inv_out_ideal'] * InverterConfig['small']['n_inv']
        
    else:
        p_ac_inv_out_preclipping = p_dc_inv_in * inv_euro_eff/100
        p_ac_inv_out_preclipping_ideal = p_dc_inv_in * inv_euro_eff/100
 
    
    RESULTS['Inverter efficiency loss'] = np.nansum(p_ac_inv_out_preclipping)/np.nansum(p_dc_inv_in)-1  # inverter loss due to efficiency
    
    
    # 22.3 CALCULATE INVERTER THRESHOLDLOSS (CLIPPING LOSS
    # calculate max inverter power output as function of temperature for actual inverter. We asusme inverter temperature is equal to air temperature (outdoor installation)
    if INVERTER['OND_import'] :
        inv_temp_ref = [T_inv_p_max, T_inv_p_nom, T_eff_1, T_eff_2]
        inv_p_ref_t_corrected = [p_ac_inv_max, p_ac_inv_nom, p_out_t1, p_out_t2]
        fill_value = (p_ac_inv_max,p_out_t2)
    
    else:
        inv_temp_ref = [T_inv_p_max, T_inv_p_nom,]
        inv_p_ref_t_corrected = [p_ac_inv_max, p_ac_inv_nom]
        fill_value = (p_ac_inv_max,p_ac_inv_nom)
    
    # create power=f(temperature) interpolation function
    inv_max_p_out_t_function = interp1d(inv_temp_ref, inv_p_ref_t_corrected, kind='slinear', bounds_error = False, fill_value=fill_value)
    
    # apply the interpolation functions above to the air temperature to find the max inverter power
    inv_max_p_out_t_corrected = inv_max_p_out_t_function(temp_air)
    p_ac_inv_out = pd.DataFrame(np.where(p_ac_inv_out_preclipping > inv_max_p_out_t_corrected*n_inv_plant, inv_max_p_out_t_corrected*n_inv_plant, p_ac_inv_out_preclipping), index = local_time)
    p_ac_inv_out_ideal = pd.DataFrame(np.where(p_ac_inv_out_preclipping_ideal > inv_max_p_out_t_corrected*n_inv_plant, inv_max_p_out_t_corrected*n_inv_plant, p_ac_inv_out_preclipping_ideal), index = local_time)
    
    # the clipping is measured on the REAL chain alone, so the row is always negative or zero
    RESULTS['Inverter clipping loss'] = np.nansum(p_ac_inv_out)/np.nansum(p_ac_inv_out_preclipping)-1
    
    # np.where above returns a DataFrame with a single unnamed column: it is turned back into a Series so that
    # the whole AC chain below stays one dimensional
    p_ac_inv_out = pd.Series(np.asarray(p_ac_inv_out).ravel(), index = local_time)
    
    
    
    # =================================================================================================
    # 22.4 - CURTAILMENT AT INVERTER OUTPUT
    # =================================================================================================
    # first of the two possible curtailment points. The AC power of the whole plant is limited right at the
    # inverter output, before the LV AC cable
    
    p_ac_plant = p_ac_inv_out
    
    if grid_limitation_point == 'inverter' :
        p_ac_plant = p_ac_inv_out.clip(upper = grid_curtail_power)
        RESULTS['Curtailment @ inverter'] = np.nansum(p_ac_plant)/np.nansum(p_ac_inv_out) - 1
        
        
    
    # =================================================================================================
    # 23 - AC-LV OHMIC LOSS
    # =================================================================================================
    
    # compute losses in cabling from inverter to trafo
    p_in_ac_lv_cable = p_ac_plant/n_inv_plant
    p_ac_lv_loss = ac_ohm_loss (p_ac = p_in_ac_lv_cable,
                                v_p2p = v_out_inv,
                                leng = ac_lv_cable_len,
                                material = 'cu',
                                section = ac_lv_cable_sect,
                                )
    
    p_out_ac_lv_cable = p_in_ac_lv_cable * (1-p_ac_lv_loss)
    RESULTS['LV AC cable loss'] = np.nansum(p_out_ac_lv_cable) / np.nansum(p_in_ac_lv_cable) - 1
    
    
    
    # =================================================================================================
    # 24 - LV/MV TRAFO LOSS
    # =================================================================================================
    
    # computes LV/MV trafo loss
    n_trafo = np.ceil(p_ac_connection/p_trafo_max)  # number of LV/MV trafo in project
    p_trafo = p_ac_connection/n_trafo               # actual trafo power
    
    p_trafo_in = p_out_ac_lv_cable * n_inv_plant / n_trafo
    
    p_trafo_out = transformer.simple_efficiency(input_power = p_trafo_in,
                                                no_load_loss = trafo_iron_loss/100,
                                                load_loss = trafo_copper_loss/100,
                                                transformer_rating = p_trafo,
                                                )
    
    RESULTS['LV/MV trafo loss'] = np.nansum(p_trafo_out)/np.nansum(p_trafo_in) - 1


    
    # =================================================================================================
    # 25 - AC-MV OHMIC LOSS
    # =================================================================================================
    
    # computes losses in MV cabling connecting trafo antennas to PV collection cabin
    
    # the auxiliaries are supplied from the MV side of the LV/MV transformer, hence their consumption is
    # subtracted between the transformer and the MV cable
    p_in_ac_mv_cable = p_trafo_out * (1+RESULTS['Auxiliaries loss'])
    
    p_ac_mv_loss = ac_ohm_loss (p_ac = p_in_ac_mv_cable,
                                v_p2p = v_connection,
                                leng = ac_mv_cable_len,
                                material = 'alu',
                                section = ac_mv_cable_sect,
                                )
    
    p_out_ac_mv_cable = p_in_ac_mv_cable * (1-p_ac_mv_loss)
    RESULTS['MV AC cable loss'] = np.nansum(p_out_ac_mv_cable) / np.nansum(p_in_ac_mv_cable) - 1
    
    
    
    # =================================================================================================
    # 26 - CONNECTION CABLE LOSS 
    # =================================================================================================
    
    # the first year degradation is applied to the hourly power BEFORE the curtailment, so that the
    # curtailment loss is computed on the power the plant will actually deliver during its first year
    p_ac_mv_cabin = p_out_ac_mv_cable * n_trafo * (1+RESULTS['First year degradation'])
    
    # CURTAILMENT AT THE ENTRY OF THE MV CONNECTION LINE
    # the plant is limited at its collection cabin, before the MV line that leaves it
    if grid_limitation_point == 'mv line in' :
        p_ac_mv_cabin_curtailed = p_ac_mv_cabin.clip(upper = grid_curtail_power)
        RESULTS['Curtailment @ MV line in'] = np.nansum(p_ac_mv_cabin_curtailed)/np.nansum(p_ac_mv_cabin) - 1
        p_ac_mv_cabin = p_ac_mv_cabin_curtailed
    
    # END OF BLOCK 1: the energy leaving the PV plant at its MV cabin
    energy_block_1 = np.nansum(p_ac_mv_cabin)
    
    # PV AVAILABILITY
    # it opens BLOCK 2 and is applied to the hourly power, so that everything downstream of it, the MV line
    # and the whole HV section, sees the power the plant actually delivers
    p_ac_mv_cabin = p_ac_mv_cabin * (1+RESULTS['PV availability'])
    
    # compute loss in connection cable leaving the PV collection parallel cabin
    p_ac_connection_cable_in = p_ac_mv_cabin / connection_cable_nr
    p_ac_connection_loss = ac_ohm_loss (p_ac = p_ac_connection_cable_in,
                                        v_p2p = v_connection,
                                        leng = connection_cable_len,
                                        material = 'alu',
                                        section = connection_cable_sect,
                                        )
    p_ac_connection_cable_out = p_ac_connection_cable_in * (1-p_ac_connection_loss)
    RESULTS['MV connection loss'] = np.nansum(p_ac_connection_cable_out) / np.nansum(p_ac_connection_cable_in) - 1
    
    p_ac_poc = p_ac_connection_cable_out * connection_cable_nr
    
    # END OF BLOCK 2: the energy at the exit of the MV connection line
    energy_block_2 = np.nansum(p_ac_poc)
    
    
    
    # =================================================================================================
    # 27 - MV/HV TRANSFORMER AND HV CONNECTION LINE
    # =================================================================================================
    # Both are optional. When no MV/HV transformer is declared, the Point of Connection sits on the MV network
    # and this whole section leaves the power untouched
    
    # CURTAILMENT AT THE ENTRY OF THE HV TRANSFORMER
    # this is the same physical node as the exit of the MV connection line
    if grid_limitation_point == 'hv trafo in' :
        p_ac_hv_trafo_in = p_ac_poc.clip(upper = grid_curtail_power)
        RESULTS['Curtailment @ HV trafo in'] = np.nansum(p_ac_hv_trafo_in)/np.nansum(p_ac_poc) - 1
        p_ac_poc = p_ac_hv_trafo_in
    
    # MV/HV TRANSFORMER LOSS
    if hv_trafo_p_nom is not None :
        p_ac_hv_trafo_out = transformer.simple_efficiency(input_power = p_ac_poc,
                                                          no_load_loss = hv_trafo_iron_loss/100,
                                                          load_loss = hv_trafo_copper_loss/100,
                                                          transformer_rating = hv_trafo_p_nom*1000*1000,
                                                          )
        
        RESULTS['HV trafo loss'] = np.nansum(p_ac_hv_trafo_out)/np.nansum(p_ac_poc) - 1
        p_ac_poc = p_ac_hv_trafo_out
    
    # CURTAILMENT AT THE EXIT OF THE HV TRANSFORMER
    # this is the same physical node as the entry of the HV line
    if grid_limitation_point == 'hv trafo out' :
        p_ac_hv_line_in = p_ac_poc.clip(upper = grid_curtail_power)
        RESULTS['Curtailment @ HV trafo out'] = np.nansum(p_ac_hv_line_in)/np.nansum(p_ac_poc) - 1
        p_ac_poc = p_ac_hv_line_in
    
    # HV CONNECTION LINE LOSS.
    # It is a flat percentage per km of line, as it is customary at this stage of a project, when the
    # conductor of the HV line is not known yet
    if hv_connection_cable_len is not None :
        hv_line_loss = hv_connection_cable_len * hv_line_loss_per_km/100
        RESULTS['HV connection loss'] = -hv_line_loss
        p_ac_poc = p_ac_poc * (1-hv_line_loss)
    
    
    
    # =================================================================================================
    # 27.1 - GRID CONNECTION LIMITATION AT THE POINT OF CONNECTION
    # =================================================================================================
    # last of the curtailment points, at the very exit of the plant. It is the useful one when a power plant
    # controller limits the whole plant instead of each single inverter
    
    if grid_limitation_point == 'poc' :
        p_ac_poc_curtailed = p_ac_poc.clip(upper = grid_curtail_power)
        RESULTS['Curtailment @ PoC'] = np.nansum(p_ac_poc_curtailed)/np.nansum(p_ac_poc) - 1
        p_ac_poc = p_ac_poc_curtailed
    
    # GRID AVAILABILITY
    # it closes BLOCK 3 and is applied to the hourly power, so that the series returned by the function is
    # the energy the plant actually injects into the grid
    p_ac_poc = p_ac_poc * (1+RESULTS['Grid availability'])
    
    # END OF BLOCK 3: the energy injected into the grid at the Point of Connection
    energy_block_3 = np.nansum(p_ac_poc)
    
    # hourly net energy delivered at the Point of Connection, in kWh. The simulation timestep is one hour,
    # hence the average power over the step is also the energy of the step
    energy_hourly = p_ac_poc/1000
    energy_hourly.name = 'Net energy at PoC [kWh]'
    
    
    
    # =================================================================================================
    # 28 - COMPUTES NET YIELD AND PERFORMANCE RATIO
    # =================================================================================================
    # The Performance Ratio is computed the same way PVsyst does, as the ratio between the specific energy
    # actually delivered and the irradiation received by the plane of array:
    #
    #     Net Yield [kWh/kWp] = energy / peak DC power
    #     PR                  = Net Yield / POAirradiation
    #
    # where POAirradiation is the global irradiation on the collector plane, front side only, before the
    # horizon and the near shading losses. That is exactly the reference PVsyst uses, which makes the two
    # Performance Ratios comparable one to one.
    #
    # NOTE: the PR is NOT the product of the loss rows. The rows describe where the energy goes, but two
    # coefficients convert the optical chain into the electrical one, the STC efficiency and the bifaciality
    # factor, and neither of them is a loss. PVsyst applies both at the 'PV conversion' step of its loss
    # diagram and computes its PR separately as Yf/Yr, which is what is reproduced here
    
    # POAirradiation on the collector plane, front side, single shed. It is the 'GlobInc' of a PVsyst report
    glob_inc = np.nansum(irradPOA_single_shed['poa_global'])/1000   # [kWh/m2]
    
    for pr_row, yield_row, block_energy in zip(pr_rows, yield_rows, [energy_block_1, energy_block_2, energy_block_3]):
        RESULTS[yield_row] = block_energy/p_dc_peak          # [Wh/Wp] = [kWh/kWp]
        RESULTS[pr_row] = RESULTS[yield_row]/glob_inc*100
    
    
    # CONSISTENCY CHECK
    # The Performance Ratio is a ratio of energies and the loss rows are only descriptive, so the two are no
    # longer forced to agree. They are still tied by an exact relation, which is worth verifying:
    #
    #     PR = product of the loss rows x bifaciality factor x inverter threshold correction
    #
    # The two corrections are the stages of the calculation that are NOT a loss row:
    #
    #  - the bifaciality factor is the optical to electrical conversion of the rear side, applied together
    #    with the STC efficiency and therefore outside the chain, exactly as PVsyst does;
    #
    # A mismatch means that a stage of the chain does not start where the previous one ended, which is exactly
    # the kind of error that stays invisible in the results table while quietly moving the numbers.
    #
    # NOTE: the bifaciality is now the ONLY correction needed. It used to take a second one, because the old
    # 'Inverter thresholds loss' row mixed the ideal and the real chain; splitting it into 'Inverter MPPT V/I loss' and
    # 'Inverter clipping loss' removed that discontinuity at the source
    bifaciality_factor = np.nansum(irradPOA_effective['poa_global'])/np.nansum(irradPOA_soiling['poa_front'] + irradPOA_soiling['poa_back'])
    
    chain_product = (RESULTS[loss_rows_all]+1).prod() * bifaciality_factor
    
    if RESULTS[pr_rows[-1]] and abs(chain_product/(RESULTS[pr_rows[-1]]/100) - 1) > 0.005:
        PVsim_err_msg += 'WARNING! The Performance Ratio computed from the delivered energy differs by more than 0.5% from the product of the loss rows corrected by the bifaciality factor. This points to a stage of the loss chain that does not start where the previous one ends\n\n'
    
    
    
    # =================================================================================================
    # 29 - CONVERT RESULTS IN % AND PREPARE OUTPUT DICTIONARY  
    # =================================================================================================
    
    #  converts results value fom unitless to %.
    RESULTS[['POA gain']] = RESULTS[['POA gain']]*100
    RESULTS[loss_rows_all] = RESULTS[loss_rows_all]*100
    RESULTS = pd.to_numeric(RESULTS, errors='coerce')
    RESULTS = RESULTS.round(2)
    
    #  concat in a DataFrame the results numerical values and the relative units
    RESULTS = pd.concat([RESULTS,results_units],axis=1)
    RESULTS.columns = ['Value','Unit']
    
    #  print result DataFrame.
    #  The options below disable every truncation: when a DataFrame holds more rows or columns than the
    #  display limits, pandas replaces the ones in the middle with '...'. option_context restores the previous
    #  settings on exit, so the user session is left untouched
    with pd.option_context('display.float_format', '{:+,.2f}'.format,
                           'display.max_rows', None,
                           'display.max_columns', None,
                           'display.width', 1000,
                           'display.expand_frame_repr', False,
                           ):
        print(RESULTS)
    

    return RESULTS, energy_hourly, PVsim_err_msg  
    
    
