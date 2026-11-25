# -*- coding: utf-8 -*-
"""
Created on Sun Oct 12 10:51:10 2025

@author: AntonioVALLONE
"""

import os
import math
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from datetime import datetime

from pvlib import tools
from pvlib import irradiance
from pvlib import bifacial

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

current_dir = os.path.dirname(os.path.abspath(__file__)) # get current file directory

def input_check_classic(POWER, LOCATION, METEO, MMS, GROUND, MODULE, INVERTER, TRAFO, LOSS, CONNECTION, CURTAILMENT):
    """
    This function does several things:
    1) Check that input for PVsim code are free of errors
    2) Infer missing/optional data
    3) if PAN or OND is provided, import the correct variables into the MODULE and INVERTER dict

    """
    errorflag = False
    check_err_msg = '\n'
    
    # CHECK POWER
    try:        
        if not(POWER['p_dc_peak']>=0):
            check_err_msg += 'ERROR! PV project DC power - invalid or missing data. Please insert valid project DC power \n\n'
            errorflag = True
    except:
            check_err_msg += 'ERROR! PV project DC power - invalid or missing data. Please insert valid project DC power \n\n'
            errorflag = True
  
    try:        
        if not(POWER['p_ac_connection']>=0):
            check_err_msg += 'ERROR! AC connection power - invalid or missing data. Please insert valid AC connection power \n\n'
            errorflag = True
    except:
            check_err_msg += 'ERROR! AC connection power - invalid or missing data. Please insert valid AC connection power \n\n'
            errorflag = True
    
    
    # CHECK LOCATION 
    try:
        if not (-90<= LOCATION['latitude'] <=90 ):
            check_err_msg += 'ERROR: Latitude - invalid or missing data. Latitude should be between -90 and 90\n\n'
            errorflag = True
    except:
        check_err_msg += 'ERROR! Latitude - invalid or missing data. Expected formats: 37.778\n\n'
        errorflag = True
    
    try:
        if not (-180<= LOCATION['longitude'] <=180 ):
            check_err_msg += 'ERROR! Longitude - invalid or missing data. Longitude should be between -180 and 180\n\n'
            errorflag = True
    except :
        check_err_msg += 'ERROR! Longitude - invalid or missing data. Expected formats: 37.778\n\n'
        errorflag = True
            
    # CHECK METEO
    try:
        if not (METEO['TMY_source'] in ['PVGIS v5.3', 'SolarGIS TMY', 'Solcast API', 'Meteonorm API']) :
            check_err_msg +='ERROR! Invalid TMY source. Valid TMY source are either: "PVGIS v5.3", "Solcast API", "SolarGIS TMY"\n\n' 
            errorflag = True
    except:
        check_err_msg += 'ERROR! Invalid TMY source. Valid TMY source are either: "PVGIS v5.3", "Solcast API", "SolarGIS TMY"\n\n'
        errorflag = True
    
    try:
        if (METEO['TMY_source'] in ['Solcast API']):
            if not (METEO['api_key'] and isinstance(METEO['api_key'], str)):
                check_err_msg += 'ERROR! If you wish to import TMY from Solcast, you must provide a valid api key\n\n' 
                errorflag = True
    except:
        check_err_msg += 'ERROR! If you wish to import TMY from Solcast, you must provide a valid api key\n\n'
        errorflag = True
        
    try:
        if (METEO['TMY_source'] in ['Meteonorm API']):
            if not (METEO['api_key'] and isinstance(METEO['api_key'], str)):
                check_err_msg += 'ERROR! If you wish to import TMY from Meteonorm, you must provide a valid api key\n\n' 
                errorflag = True
    except:
        check_err_msg += 'ERROR! If you wish to import TMY from Solcast, you must provide a valid api key\n\n'
        errorflag = True
        
    try:
        if (METEO['TMY_source'] in ['SolarGIS TMY']):
            if not (METEO['TMY_file_name'] and isinstance(METEO['TMY_file_name'], str)):
                check_err_msg += 'ERROR! If you wish to import TMY from your own SolarGIS TMY, you must provide a valid TMY location\n\n'
                errorflag = True
    except:
        check_err_msg += 'ERROR! Invalid TMY source. Valid TMY source are either: "PVGIS v5.3", "Solcast API", "SolarGIS TMY"\n\n'
        errorflag = True
        
    # CHECK GROUND
    try:
        if not (0<=GROUND['albedo']<=1):
            check_err_msg += 'WARNING! Albedo - invalid or missing data. Albedo must be between 0 and 1. The simulation assumed an albedo of 0.2 (=20%)\n\n'
            GROUND['albedo'] = 0.2
    except:
        check_err_msg += 'WARNING! Albedo - invalid or missing data. Albedo must be between 0 and 1. The simulation assumed an albedo of 0.2 (=20%)\n\n'
        GROUND['albedo'] = 0.2

    try:
        if not (0<=GROUND['slope_tilt']<=1000):
            check_err_msg += 'WARNING! Ground slope - invalid or missing data. Ground slope should be expressed in %, from 0% (flat) to 1000% (vertical wall). The simulation assumed a flat ground (ground slope = 0%)\n\n'
            GROUND['slope_tilt'] = 0
    except:
        check_err_msg += 'WARNING! Ground slope - invalid or missing data. Ground slope should be expressed in degrees %, from 0% (flat) to 1000% (vertical wall). The simulation assumed a flat ground (ground slope = 0%)\n\n'
        GROUND['slope_tilt'] = 0
    
    try:
        if not (0<=GROUND['slope_azimuth']<=360):
            check_err_msg += 'WARNING! Ground slope azimuth format - invalid or missing data. Ground slope azimuth is the direction where a ball would roll if left on the slope (0 = North, 90 = East, 180 = South, 270 = West). It must be expressed in degrees °. The simulation assumed a flat ground\n\n'
            GROUND['slope_tilt'] = 0
            GROUND['slope_azimuth'] = 0
    except:
        check_err_msg += 'WARNING! Ground slope azimuth format - invalid or missing data. Ground slope azimuth is the direction where a ball would roll if left on the slope (0 = North, 90 = East, 180 = South, 270 = West). It must be expressed in degrees °. The simulation assumed a flat ground\n\n'
        GROUND['slope_tilt'] = 0
        GROUND['slope_azimuth'] = 0
        
        
    # CHECK PV MODULES
    # extract data that will be lost with PAN import. Such data are not present in the PAN file, therefore will be verwritten and lost if the PAN data is imported
    PAN_import = MODULE['PAN_import']
    
    
    # Import module data from PAN
    panerrorflag = False
    
    try:
        if not (isinstance(MODULE['PAN_import'], bool)):
            check_err_msg += 'ERROR! Not clear if I need to import a PV module data from a PAN. Check PAN import option\n\n'
            errorflag = True
            panerrorflag = True
    except:
        check_err_msg += 'ERROR! Not clear if I need to import a PV module data from a PAN. Check PAN import option\n\n'
        errorflag = True
        panerrorflag = True
    
    try:  
        if MODULE['PAN_import']:
            if not (MODULE['PAN_file_name'] and isinstance(MODULE['PAN_file_name'], str)):
                check_err_msg += 'ERROR! If you wish to import PV module data from PAN, you must provide a valid PAN file\n\n'
                errorflag = True
                panerrorflag = True
            else:
                PAN_path = os.path.join(current_dir, 'PAN', MODULE['PAN_file_name'])
                MODULE = pan_import_classic(PAN_path)
                
                #restore MODULE data deleted during import of OND
                MODULE['PAN_import'] = PAN_import

    except:
        check_err_msg += 'ERROR! Could not import PV module data from PAN provided\n\n'
        errorflag = True
        panerrorflag = True
    
    # CHECK MODULE DATA
    if not panerrorflag:            
        
        try:
            float(MODULE['PNom'])
        except: 
            check_err_msg += 'ERROR! PV module power - invalid or missing data\n\n'
            errorflag = True
            
        try:
            float(MODULE['Height'])
        except: 
            check_err_msg += 'ERROR! PV module height - invalid or missing data\n\n'
            errorflag = True
    
        try:
            float(MODULE['Width'])
        except: 
            check_err_msg += 'ERROR! PV module width - invalid or missing data\n\n'
            errorflag = True    
    
        try:
            if not(isinstance(MODULE['NCelS'], int)):
                check_err_msg += 'WARNING! PV module cells in series must be an integer. The simulation assumed a default value of 66\n\n'
                MODULE['NCelS'] = 66
        except: 
            check_err_msg += 'WARNING! PV module cells in series must be an integer. The simulation assumed a default value of 66\n\n'
            MODULE['NCelS'] = 66
    
        try:
            if not(isinstance(MODULE['NCelP'], int)):
                check_err_msg += 'WARNING! PV module cells in parallel must be an integer. The simulation assumed a default value of 2\n\n'
                MODULE['NCelP'] = 2
        except: 
            check_err_msg += 'WARNING! PV module cells in parallel must be an integer. The simulation assumed a default value of 2\n\n'
            MODULE['NCelP'] = 2
    
        try:
            if not(isinstance(MODULE['NDiode'], int)):
                check_err_msg += 'WARNING! Number of bypass diodes in PV module is invalid. The simulation assumed a default value of 3\n\n'
                MODULE['NDiode'] = 3
        except: 
            check_err_msg += 'WARNING! Number of bypass diodes in PV module is invalid. The simulation assumed a default value of 3\n\n'
            MODULE['NDiode'] = 3
    
        try:
            float(MODULE['Imp'])
        except: 
            check_err_msg += 'ERROR! PV module mppt current - invalid or missing data\n\n'
            errorflag = True
            
        try:
            float(MODULE['Voc'])
        except: 
            check_err_msg += 'ERROR! PV module open circuit voltage - invalid or missing data\n\n'
            errorflag = True    
            
        try:
            float(MODULE['Vmp'])
        except:
            check_err_msg += 'WARNING! PV module mppt voltage - invalid or missing data. Software assumed a default value as Pnom/Imp\n\n'
            MODULE['Vmp'] = MODULE['Pnom']/MODULE['Imp']
                
        try:
            float(MODULE['muISC'])
        except:
            check_err_msg += 'WARNING! PV module short circuit temperature coefficient - invalid or missing data. The simulation assumed a default value of 6.75\n\n'
            MODULE['muISC'] = 6.75
    
        try:
            float(MODULE['Gamma'])
        except:
            check_err_msg += 'WARNING: PV module Gamma (Ideality factor) - invalid or missing data. The simulation assumed a default value of 1.02\n\n'
            MODULE['Gamma'] = 1.02
    
        try:
            float(MODULE['muGamma'])   
        except:
            check_err_msg += 'WARNING! PV module muGamma (temperature coefficient for ideality factor ) - invalid or missing data. The simulation assumed a default value of -0.0002\n\n'
            MODULE['muGamma'] = -0.0002
    
        try:
            float(MODULE['Isc'])   
        except:
            check_err_msg += 'WARNING! PV module short circuit current - invalid or missing data. The simulation assumed a value derived from other MODULE data\n\n'
            MODULE['Isc'] = MODULE['Imp'] + 1
    
        try:
            float(MODULE['RShunt'])
        except:
            check_err_msg += 'WARNING! PV module shunt resistance - invalid or missing data. The simulation assumed a value derived from other MODULE data\n\n'
            MODULE['RShunt'] = MODULE['Vmp']/(0.2 * (MODULE['Isc'] - MODULE['Imp']))   # from PVsyst guide
    
        try:
            float(MODULE['Rp_0'])
        except:
            check_err_msg += 'WARNING! PV module shunt resistance at null irradiance (Rp_0) - invalid or missing data. The simulation assumed a default value of 4xRShunt\n\n'
            MODULE['Rp_0'] = MODULE['RShunt'] * 4     # from PVsyst guide
            
        try:
            float(MODULE['Rp_Exp'])
        except:
            check_err_msg += 'WARNING! PV module exponential factor (Rp_Exp) - invalid or missing data. The simulation assumed a default value of 5.5\n\n'
            MODULE['Rp_Exp'] = 5.5
            
        try:
            float(MODULE['RSerie'])
        except:
            check_err_msg += 'WARNING! PV module serie resistance exponent (RSerie) - invalid or missing data. The simulation will try to assume a value of derived from other MODULE data\n\n'
            MODULE['RSerie'] = 8/(MODULE['Isc']*3)      # made up formula
            
            """
            About RSerie calculation
            Rserie is the slope of the IV curve at Voc. In an ideal diode is 0
            It can be overestimated as (Voc-Vmp)/(Imp) (slope between Voc and Max power point)
            
            To match the slope, I found out that this formula works better: (Voc-Vmp)/(3xImp)
            If Voc is not provided, then we assume that Voc-Vmp = 8 V 
            
            """
        
        try:
            float(MODULE['Absorb'])
        except:
            check_err_msg += 'WARNING! PV module absorption coefficient - invalid or missing data. The simulation assumed a default value of 0.90\n\n'
            MODULE['Absorb'] = 0.9     # from PVsyst guide       
    
        try:
            float(MODULE['BifacialityFactor'])
        except:
            check_err_msg += 'WARNING! Bifaciality factor - invalid or missing data. The simulation assumed a default value of 0.80\n\n'
            MODULE['BifacialityFactor'] = 0.8      # from PVsyst guide 
    
        try:
            if not ('half' in str.lower(MODULE['SubModuleLayout']) or 'InLength' in str.lower(MODULE['SubModuleLayout'])):
                MODULE['SubModuleLayout'] = 'slTwinHalfCells'
                check_err_msg += "ERROR! Invalid PV module substring layout. Module layout should be 'TwinHalfCells' (Half cut cells) or 'FullCells' (full cut cells). The simulation assumed halfcut module\n\n"
        except:
            MODULE['SubModuleLayout'] = 'TwinHalfCells'
            check_err_msg += "ERROR! Invalid PV module substring layout. Module layout should be 'TwinHalfCells' (Half cut cells) or 'FullCells' (full cut cells). The simulation assumed halfcut module\n\n"

    if LOSS['iam_model'] =='PAN':
        try:
            if MODULE['iam_theta'] is None:
                MODULE['iam_theta'] = [0, 30, 50, 60, 70, 75, 80, 85, 90]
                check_err_msg += "WARNING! Invalid PV iam reference angle data. The simulation assumed default iam reference angle data\n\n"
        except: 
            MODULE['iam_theta'] = [0, 30, 50, 60, 70, 75, 80, 85, 90]
            check_err_msg += "WARNING! Invalid PV iam reference angle data. The simulation assumed default iam reference angle data\n\n"
        
        try:
            if MODULE['iam_loss'] is None:
                MODULE['iam_loss'] = [1, 1, 0.999, 0.994, 0.969, 0.928, 0.829, 0.585, 0]
                check_err_msg += "WARNING! Invalid PV iam reference loss data. The simulation assumed default iam reference loss data\n\n"
        except: 
            MODULE['iam_loss'] = [1, 1, 0.999, 0.994, 0.969, 0.928, 0.829, 0.585, 0]
            check_err_msg += "WARNING! Invalid PV iam reference loss data. The simulation assumed default iam reference loss data\n\n"
        
        if len(MODULE['iam_loss']) != len(MODULE['iam_theta']):
            check_err_msg += "ERROR! Invalid PV iam data. iam reference angle and iam refence loss must have same number of elements\n\n"
            errorflag = True
    
    try:
        float(MODULE['VRevDiode'])
    except:
        check_err_msg += 'WARNING! PV module bypass diode voltage drop - invalid or missing data. The simulation assumed a default value of 0.7\n\n'
        MODULE['VRevDiode'] = 0.7     

# CHECK MODULE MOUNTING STRUCTURE
    try:
        if MMS['mount_type'] == 'Single-axis tracker':
            angle = 0
        elif MMS['mount_type'] == 'Fixed tilt':
            angle = MMS['tilt_angle']
        else:
            check_err_msg += 'ERROR! Invalid mouting mount type. Mount type can be either "Single-axis tracker" (Single Axis Tracker) or "Fixed tilt" (Fixed Tilt Structure)\n\n'
            errorflag = True
    except:
        check_err_msg += 'ERROR! Invalid mouting mount type. Mount type can be either "Single-axis tracker" (Single Axis Tracker) or "Fixed tilt" (Fixed Tilt Structure)\n\n'
        errorflag = True
    
    try:
        if not ((MMS['mod_configuration'] in ['1P', '2P', '3P','4P','1L','2L','3L','4L','5L','6L'])):
            check_err_msg += "ERROR! Invalid module configuration.Accepted configurations are: '1P', '2P', '3P', '4P' (Portrait) or '1L', '2L', '3L', '4L', '5L', '6L' (Landscape)\n\n"
            errorflag = True
    except:
        check_err_msg += "ERROR! Invalid module configuration. Accepted configurations are: '1P', '2P', '3P', '4P' (Portrait) or '1L', '2L', '3L', '4L', '5L', '6L' (Landscape)\n\n"
        errorflag = True
    
    try:
        if not isinstance(MMS['modules_per_string'], int):
            check_err_msg += 'ERROR! invalid number of modules per string. This must be an integer\n\n'
            errorflag = True
    except:
        check_err_msg += 'ERROR! Modules per string - invalid or missing data\n\n'
        errorflag = True
    
    if not errorflag:
        if MMS['mod_configuration'][1] =='P':
            min_pitch = MODULE['Height'] * int(MMS['mod_configuration'][0]) * tools.cosd(angle)
        elif MODULE['mod_configuration'][1] =='L':
            min_pitch = MODULE['Width'] * int(MMS['mod_configuration'][0]) * tools.cosd(angle)
    else: min_pitch = 0

    try:    
        if MMS['pitch'] < min_pitch:
            check_err_msg += 'ERROR! the pitch you selected is incompatible with the configuration and module size. The pitch should be at least ' + str(int(math.ceil(min_pitch*10)/10)) + 'm\n\n'
            errorflag = True
    except: 
        check_err_msg += 'ERROR! Pitch - invalid or missing data\n\n'
        errorflag = True
    
    try:
        if not (0<=MMS['tilt_angle']<=90):
            check_err_msg += 'ERROR! Invalid tilt angle - invalid or missing data. Tilt angle should be between 0 and 90 degrees\n\n'
            errorflag = True
    except:
        check_err_msg += 'ERROR! Invalid tilt angle - invalid or missing data. Tilt angle should be between 0 and 90 degrees\n\n'
        errorflag = True

    try:
        if not (0<=MMS['axis_azimuth']<=360):
            if MMS['mount_type'] == 'Single-axis tracker':
                check_err_msg += 'WARNING! Structure axis azimuth - invalid or missing data. Axis azimuth angle should be between 0 and 360 degrees. The simulation assumed a value of 0°, compatible with Tracker structure oriented along N-S axis\n\n'
                MMS['axis_azimuth'] = 180
            elif MMS['mount_type'] == 'Fixed tilt':
                check_err_msg += 'WARNING! Structure axis azimuth - invalid or missing data. Axis azimuth angle should be between 0 and 360 degrees. The simulation assumed a value of 90°, compatible with fixed tilt structure oriented along E-W axis\n\n'
                MMS['axis_azimuth'] = 90
    except:
        try:
            if MMS['axis_azimuth']== None or MMS['axis_azimuth']=='':
                if MMS['mount_type'] == 'Single-axis tracker':
                    check_err_msg += 'WARNING! Structure axis azimuth - invalid or missing data. Axis azimuth angle should be between 0 and 360 degrees. The simulation assumed a value of 0°, compatible with Tracker structure oriented along N-S axis\n\n'
                    MMS['axis_azimuth'] =180
                elif MMS['mount_type'] == 'Fixed tilt':
                    check_err_msg += 'WARNING! Structure axis azimuth - invalid or missing data. Axis azimuth angle should be between 0 and 360 degrees. The simulation assumed a value of 90°, compatible with fixed tilt structure oriented along E-W axis\n\n'
                    MMS['axis_azimuth'] =90
        except:
            check_err_msg += 'ERROR! Structure axis azimuth - invalid or missing data. Axis azimuth angle should be between 0 and 360 degrees\n\n'
            errorflag = True

    try:
        if not isinstance(MMS['backtrack'], bool):
            check_err_msg += 'ERROR! Invalid backtracking format\n\n'
            errorflag = True
    except:
            check_err_msg += 'ERROR! Invalid backtracking format\n\n'
            errorflag = True
    
    try:
        float(MMS['ground_gap'])
    except:
        check_err_msg += 'WARNING! ground clearance - invalid or missing data. The simulation assumed a value of 0.5m\n\n'
        MMS['ground_gap'] = 0.5

    # INVERTER CHECK
    # extract data that will be lost with OND import. Such data are not present in the OND file, therefore will be rewritten and lost if the OND data is imported
    OND_import = INVERTER['OND_import']
    
    # Import inverter data from OND
    onderrorflag = False
    
    try:
        if not (isinstance(INVERTER['OND_import'], bool)):
            check_err_msg += 'ERROR! Not clear if I need to import inverter data from a OND. Check OND import option\n\n'
            errorflag = True
            onderrorflag = True
    except:
        check_err_msg += 'ERROR! Not clear if I need to import inverter data from a OND. Check OND import option\n\n'
        errorflag = True
        onderrorflag = True
    
    try:  
        if INVERTER['OND_import']:
            if not (INVERTER['OND_file_name'] and isinstance(INVERTER['OND_file_name'], str)):
                check_err_msg += 'ERROR! If you wish to import inverter data from OND, you must provide a valid OND file\n\n'
                errorflag = True
                onderrorflag = True
            else:
                OND_path = os.path.join(current_dir, 'OND', INVERTER['OND_file_name'])
                INVERTER = ond_import_classic(OND_path)
                
                #restore INVERTER data deleted during import of OND
                INVERTER['OND_import'] = OND_import

    except:
        check_err_msg += 'ERROR! Could not import inverter data from OND provided\n\n'
        errorflag = True
        onderrorflag = True
    
  # check other inverter data
    if not onderrorflag:
    
        try: 
            float(INVERTER['PMaxOUT'])
        except:
            check_err_msg += 'ERROR! inverter max AC output power - invalid or missing data\n\n'
            errorflag = True
            
        try: 
            float(INVERTER['PMaxOUT'])
        except:
            check_err_msg += 'ERROR! inverter max AC output power - invalid or missing data\n\n'
            errorflag = True
    
        try: 
            float(INVERTER['VMppMin'])
        except:
            check_err_msg += 'ERROR! inverter minimum MPPT voltage - invalid format or not specified\n\n'
            errorflag = True
            
        try: 
            float(INVERTER['VMPPMax'])
        except:
            check_err_msg += 'ERROR! Inverter maximum MPPT voltage - invalid or missing data\n\n'
            errorflag = True
    
    
        if INVERTER['OND_import']:
            try:
                INVERTER['i_max_mppt'] = INVERTER['IMaxDC']/INVERTER['NbMPPT']
            except:
                check_err_msg += 'ERROR! Inverter max current per MPPT \'IMaxDC\' or inverter number of inputs \'NbInputs\' not present on OND. Cannot compute essnetial data for simulation. Check your OND\n\n'
                errorflag = True
            
        else:
            try:
                float(INVERTER['i_max_mppt'])
            except:
                check_err_msg += 'ERROR! Inverter max current per MPPT - invalid or missing data\n\n'
                errorflag = True
    
        try:
            isinstance(INVERTER['NbMPPT'], int)
        except:
            check_err_msg += 'ERROR! Inverter total number of input string - invalid or missing data\n\n'
            errorflag = True 
    
        try:
            isinstance(INVERTER['NbInputs'], int)
        except:
            check_err_msg += 'ERROR! Inverter total number of input string - invalid or missing data\n\n'
            errorflag = True 
            
        try: 
            float(INVERTER['VOutConv'])
        except:
            check_err_msg += 'WARNING! Inverter output AC voltage - invalid or missing data. The simulation assumed a value of 800 V\n\n'
            INVERTER['VOutConv'] = 800
    
        
        if not INVERTER['OND_import']:
            try:
                if not 0<=float(INVERTER['inv_euro_eff'])<=100 :
                  check_err_msg += 'WARNING! Inverter european efficiency - invalid or missing data. Efficiency should be between 0 and 100 [%]. The simulation assumed a value of of 98.8%\n\n'
            except:
                check_err_msg += 'WARNING! Inverter european efficiency - invalid or missing data. Efficiency should be between 0 and 100 [%]. The simulation assumed a value of 98.8%\n\n'

        try: 
            float(INVERTER['TPNom'])
        except:
            check_err_msg += 'WARNING! Inverter nominal temperature - invalid or missing data. The simulation assumed a value of 40 °C\n\n'
            INVERTER['TPNom'] = 40

        try: 
            float(INVERTER['TPMax'])
        except:
            check_err_msg += 'WARNING! Inverter max power temperature - invalid or missing data. The simulation assumed a value of 30 °C\n\n'
            INVERTER['TPMax'] = 30


    # TRAFO input check
    try:
        float(TRAFO['p_trafo_max'])
    except:
        check_err_msg += 'WARNING: Max trafo power - invalid or missing data. The simulation assumed a value of 6.4MVA\n\n'
        TRAFO['p_trafo_max'] = 6.4


    # LOSSES input check
    try:
        if not (LOSS['iam_model'] in ['PAN' , 'GlassAR']) :
            check_err_msg +='WARNING! Invalid iam loss model. Valid iam model are either: "PAN" or "GlassAR". For this simulation it will be assumed "GlassAR" model \n\n' 
            LOSS['iam_model'] = 'GlassAR'
    except:
        check_err_msg +='WARNING! Invalid iam model. Valid iam loss model are either: "PAN" or GlassAR". For this simulation it will be assumed "GlassAR" model \n\n' 
        LOSS['iam_model'] = 'GlassAR'
    
    try:
        if not (0<= LOSS['soiling_loss'] <=100 ):
            check_err_msg += 'WARNING! Soiling loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 2%\n\n'
            LOSS['soiling_loss'] = 2
    except:
        check_err_msg += 'WARNING! Soiling loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 2%\n\n'
        LOSS['soiling_loss'] = 2
        
    try:
        if not (-100<= LOSS['module_quality'] <=100 ):
            check_err_msg += 'WARNING! Module quality - invalid or missing data. Loss should be between -100 and -100 [%]. The simulation assumed a value of -0.2%\n\n'
            LOSS['module_quality'] = -0.2
    except:
        check_err_msg += 'WARNING! Module quality - invalid or missing data. Loss should be between -100 and -100 [%]. The simulation assumed a value of -0.2%\n\n'
        LOSS['module_quality'] = -0.2

    try:
        if not (0<= LOSS['lid_loss'] <=100 ):
            check_err_msg += 'WARNING! LID loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 0.3%\n\n'
            LOSS['lid_loss'] = 0.3
    except:
        check_err_msg += 'WARNING! LID loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 0.3%\n\n'
        LOSS['lid_loss'] = 0.3   

    try:
        if not (0<= LOSS['mismatch_loss_front'] <=100 ):
            check_err_msg += 'WARNING! Front mismatch loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 0.6%\n\n'
            LOSS['mismatch_loss_front'] = 0.6
    except:
        check_err_msg += 'WARNING! Front mismatch loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 0.6%\n\n'
        LOSS['mismatch_loss_front'] = 0.6

    try:
        if not (0<= LOSS['mismatch_loss_back'] <=100 ):
            check_err_msg += 'WARNING! Back mismatch loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 0.2%\n\n'
            LOSS['mismatch_loss_back'] = 0.2
    except:
        check_err_msg += 'WARNING! Back mismatch loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 0.2%\n\n'
        LOSS['mismatch_loss_back'] = 0.2    

    try:
        if not (0<= LOSS['dc_ohmic_loss_STC'] <=100 ):
            check_err_msg += 'WARNING! DC ohmic loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 1.5%\n\n'
            LOSS['dc_ohmic_loss_STC'] = 1.5
    except:
        check_err_msg += 'WARNING! DC ohmic loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 1.5%\n\n'
        LOSS['dc_ohmic_loss_STC'] = 1.5

    try: 
        if not (0<= LOSS['auxiliaries_loss'] <=100 ):
            try:
                if MMS['mount_type'] =='Single-axis tracker':
                    default_aux_loss = 0.4
                elif MMS['mount_type'] =='Fixed tilt':
                    default_aux_loss = 0.2
                check_err_msg += 'WARNING! Auxiliaries loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of ' + str(default_aux_loss) + '%, compatible with ' + MMS['mount_type'] + ' mounting type\n\n'
                LOSS['auxiliaries_loss'] = default_aux_loss
            except: 
                errorflag = True  # if an error is detected, is because there is an invalid mounting type. The error message has been already delivered to user
    except:
        try:
            if MMS['mount_type'] =='Single-axis tracker':
                default_aux_loss = 0.4
            elif MMS['mount_type'] =='Fixed tilt':
                default_aux_loss = 0.2
            check_err_msg +='WARNING! Auxiliaries loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of ' + str(default_aux_loss) + '%, compatible with ' + MMS['mount_type'] + ' mounting type\n\n'
            LOSS['auxiliaries_loss'] = default_aux_loss
        except: 
            errorflag = True  # if an error is detected, is because there is an invalid mounting type. The error message has been already delivered to user

    try:
        if not (0<= LOSS['trafo_iron_loss'] <=100 ):
            check_err_msg +='WARNING! Trasfromer iron loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 0.1%\n\n'
            LOSS['trafo_iron_loss'] = 0.1
    except:
        check_err_msg +='WARNING! Trasfromer iron loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 0.1%\n\n'
        LOSS['trafo_iron_loss'] = 0.1

    try:
        if not (0<= LOSS['trafo_copper_loss'] <=100 ):
            check_err_msg +='WARNING! Trasfromer copper loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 1.0%\n\n'
            LOSS['trafo_copper_loss'] = 1.0
    except:
        check_err_msg +='WARNING! Trasfromer copper loss - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 1.0%\n\n'
        LOSS['trafo_copper_loss'] = 1.0


    try: 
        if not (0<= LOSS['pv_availability'] <=100 ):
            try:
                if MMS['mount_type'] =='Single-axis tracker':
                    default_availability = 99
                elif MMS['mount_type'] =='Fixed tilt':
                    default_availability = 99.5
                check_err_msg += 'WARNING! PV Availability - invalid or missing data. Availability should be between 0 and 100 [%]. The simulation assumed a value of ' + str(default_availability) + '%, compatible with ' + MMS['mount_type'] + ' mounting type\n\n'
                LOSS['pv_availability'] = default_availability
            except: 
                errorflag = True  # if an error is detected, is because there is an invalid mounting type. The error message has been already delivered to user
    except:
        try:
            if MMS['mount_type'] =='Single-axis tracker':
                default_availability = 99
            elif MMS['mount_type'] =='Fixed tilt':
                default_availability = 99.5
            check_err_msg += 'WARNING! PV Availability - invalid or missing data. Availability should be between 0 and 100 [%]. The simulation assumed a value of ' + str(default_availability) + '%, compatible with ' + MMS['mount_type'] + ' mounting type\n\n'
            LOSS['pv_availability'] = default_availability
        except: 
            errorflag = True  # if an error is detected, is because there is an invalid mounting type. The error message has been already delivered to user

    try:
        if not (0<= LOSS['grid_availability'] <=100 ):
            check_err_msg+= 'WARNING! Grid Availability - invalid or missing data. Availability should be between 0 and 100 [%]. The simulation assumed a value of 99.7%\n\n'
            LOSS['grid_availability'] = 99.7
    except:
        check_err_msg+= 'WARNING! Grid Availability - invalid or missing data. Availability should be between 0 and 100 [%]. The simulation assumed a value of 99.7%\n\n'
        LOSS['grid_availability'] = 99.7

    try:
        if not (0<= LOSS['first_year_degrad'] <=100 ):
            check_err_msg+= 'WARNING! Plant First year degradation - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 0.5%\n\n'
            LOSS['first_year_degrad'] = 0.5
    except:
        check_err_msg+= 'WARNING! First year degradation - invalid or missing data. Loss should be between 0 and 100 [%]. The simulation assumed a value of 0.5%\n\n'
        LOSS['first_year_degrad'] = 0.5

    try:
        if not (LOSS['thermal_loss_factor']>=0):
            check_err_msg+= 'WARNING! Thermal loss factor - invalid or missing data. The simulation assumed a value of 29 W/m2/K\n\n'
            LOSS['thermal_loss_factor'] = 29
    except: 
        check_err_msg+= 'WARNING! Thermal loss factor - invalid or missing data. The simulation assumed a value of 29 W/m2/K\n\n'
        LOSS['thermal_loss_factor'] = 29

    try:
        if not (0<=LOSS['rear_shade_factor']<=100):
            check_err_msg += 'WARNING! Rear shade factor - invalid or missing data. Value should be between 0 and 100 [%]. The simulation assumed a value of 5%\n\n'
            LOSS['rear_shade_factor'] = 5 
    except:
        check_err_msg += 'WARNING! Rear shade factor - invalid or missing data. Value should be between 0 and 100 [%]. The simulation assumed a value of 5%\n\n'
        LOSS['rear_shade_factor'] = 5         

    try:
        if not(LOSS['ac_lv_cable_len']>=0):
            check_err_msg+= 'WARNING! AC LV cable lenght (from inverter to Trafo) - invalid or missing data. The simulation assumed a value of 500 m\n\n'
            LOSS['ac_lv_cable_len'] = 500
    except: 
        check_err_msg+= 'WARNING! AC LV cable lenght (from inverter to Trafo) - invalid or missing data. The simulation assumed a value of 500 m\n\n'
        LOSS['ac_lv_cable_len'] = 500

    try:
        if not(LOSS['ac_lv_cable_sect']>=0):
            check_err_msg+= 'WARNING! AC LV cable section (from inverter to Trafo) - invalid or missing data. The simulation assumed a value of 400 mm^2, compatible with most string inverters \n\n'
            LOSS['ac_lv_cable_sect'] = 400
    except: 
        check_err_msg+= 'WARNING! AC LV cable section (from inverter to Trafo) - invalid or missing data. The simulation assumed a value of 400 mm^2, compatible with most string inverters \n\n'
        LOSS['ac_lv_cable_sect'] = 400        

    try:
        if not(LOSS['ac_mv_cable_len']>=0):
            check_err_msg+= 'WARNING! AC MV cable lenght (from trafo to collection abin) - invalid or missing data. The simulation assumed a value of 1000 m\n\n'
            LOSS['ac_mv_cable_len'] = 1000
    except: 
        check_err_msg+= 'WARNING! AC LV cable lenght (from trafo to collection cabin) - invalid or missing data. The simulation assumed a value of 1000 m\n\n'
        LOSS['ac_mv_cable_len'] = 1000

    try:
        if not(LOSS['ac_mv_cable_sect']>=0):
            check_err_msg+= 'WARNING! AC LV cable section (from trafo to collection cabin) - invalid or missing data. The simulation assumed a value of 300 mm^2 \n\n'
            LOSS['ac_mv_cable_sect'] = 400
    except: 
        check_err_msg+= 'WARNING! AC LV cable section (from trafo to collection cabin) - invalid or missing data. The simulation assumed a value of 300 mm^2 \n\n'
        LOSS['ac_mv_cable_sect'] = 400       


    # CONNECTION CHECK
    # NOTE: a trailing comma on an assignment turns the value into a TUPLE, e.g. writing
    #       CONNECTION['hv_trafo_p_nom'] = None,   makes the value (None,) and not None.
    # Every optional entry is therefore unwrapped first, so that a stray comma in the launcher does not turn
    # into an obscure failure several hundred lines later
    for connection_key in list(CONNECTION.keys()):
        if isinstance(CONNECTION[connection_key], tuple) and len(CONNECTION[connection_key]) == 1:
            check_err_msg+= 'WARNING! CONNECTION[\'' + connection_key + '\'] was given as a tuple, most likely because of a trailing comma in the launcher. The single value it contains has been used \n\n'
            CONNECTION[connection_key] = CONNECTION[connection_key][0]

    try:
        if not(CONNECTION['mv_v_connection']>=0):
             check_err_msg+= 'WARNING! MV connection voltage - invalid or missing data. The simulation assumed a value of 20kV \n\n'
             CONNECTION['mv_v_connection'] = 20
    except:
        check_err_msg+= 'WARNING! MV connection voltage - invalid or missing data. The simulation assumed a value of 20kV \n\n'
        CONNECTION['mv_v_connection'] = 20
            
                          
    try:
        if not(CONNECTION['mv_connection_cable_len']>=0):
            check_err_msg+= 'WARNING! MV connection cable lenght - invalid or missing data. The simulation assumed no MV connection cable loss \n\n'
            CONNECTION['mv_connection_cable_len'] = None
    except:
       check_err_msg+= 'WARNING! MV connection cable lenght - invalid or missing data. The simulation assumed no MV connection cable loss \n\n'
       CONNECTION['mv_connection_cable_len'] = None          

    try:
        if not(CONNECTION['mv_connection_cable_sect']>=0):
            check_err_msg+= 'WARNING! MV connection cable section - invalid or missing data. The simulation assumed a standard connection cable loss \n\n'
            CONNECTION['mv_connection_cable_sect'] = None   
    except:
        check_err_msg+= 'WARNING! MV connection cable section - invalid or missing data. The simulation assumed a standard connection cable loss \n\n'
        CONNECTION['mv_connection_cable_sect'] = None   

    try:
        if not(CONNECTION['mv_connection_cable_nr']>=0):
            check_err_msg+= 'WARNING! MV connection cable number - invalid or missing data. The simulation assumed a standard connection cable loss \n\n'
            CONNECTION['mv_connection_cable_nr'] = None   
    except:
        check_err_msg+= 'WARNING! MV connection cable number - invalid or missing data. The simulation assumed a standard connection cable loss \n\n'
        CONNECTION['mv_connection_cable_nr'] = None   


    # HV TRAFO CHECK
    # The MV/HV transformer is optional: when its rated power is missing, the Point of Connection sits on the
    # MV network and both the transformer and the HV line are skipped
    try:
        if not(CONNECTION['hv_trafo_p_nom']>0):
            CONNECTION['hv_trafo_p_nom'] = None
    except:
        CONNECTION['hv_trafo_p_nom'] = None

    if CONNECTION['hv_trafo_p_nom'] is None:
        check_err_msg+= 'WARNING! MV/HV transformer rated power not given. The simulation assumed no MV/HV transformer, i.e. the Point of Connection sits on the MV network \n\n'
        CONNECTION['hv_trafo_iron_loss'] = None
        CONNECTION['hv_trafo_copper_loss'] = None
    else:
        if CONNECTION['hv_trafo_p_nom'] < POWER['p_ac_connection']:
            check_err_msg+= 'ERROR! The MV/HV transformer rated power is lower than the AC power of the plant \n\n'
            errorflag = True

        try:
            if not(0 <= CONNECTION['hv_trafo_iron_loss'] <= 100):
                raise ValueError
        except:
            check_err_msg+= 'WARNING! MV/HV transformer iron loss - invalid or missing data. The simulation assumed a value of 0.4% \n\n'
            CONNECTION['hv_trafo_iron_loss'] = 0.4

        try:
            if not(0 <= CONNECTION['hv_trafo_copper_loss'] <= 100):
                raise ValueError
        except:
            # on a large MV/HV transformer the copper loss expressed in % of the rating decreases with the size
            CONNECTION['hv_trafo_copper_loss'] = 0.04 if CONNECTION['hv_trafo_p_nom'] >= 20 else 0.02
            check_err_msg+= 'WARNING! MV/HV transformer copper loss - invalid or missing data. The simulation assumed a value of ' + str(CONNECTION['hv_trafo_copper_loss']) + '% \n\n'

        # on a large MV/HV transformer the iron loss is the dominant one when expressed in % of the rating,
        # which is the opposite of what happens on the small LV/MV transformers of the PV field
        if CONNECTION['hv_trafo_copper_loss'] > CONNECTION['hv_trafo_iron_loss']:
            check_err_msg+= 'WARNING! The MV/HV transformer copper loss is higher than its iron loss. On a large transformer the opposite is expected. Please check the inputs \n\n'


    # HV CONNECTION LINE CHECK
    try:
        if not(CONNECTION['hv_connection_cable_len']>=0):
            CONNECTION['hv_connection_cable_len'] = None
    except:
        CONNECTION['hv_connection_cable_len'] = None

    if CONNECTION['hv_connection_cable_len'] is None:
        check_err_msg+= 'WARNING! HV connection line lenght not given. The simulation assumed no HV connection line \n\n'
    elif CONNECTION['hv_trafo_p_nom'] is None:
        check_err_msg+= 'WARNING! An HV connection line was declared without a MV/HV transformer. The HV line has been ignored \n\n'
        CONNECTION['hv_connection_cable_len'] = None


    # CURTAILMENT CHECK
    # the value is normalised (lowercase, single spaces) so that the user input is NOT case sensitive
    valid_curtail_point = ['inverter', 'mv line in', 'hv trafo in', 'hv trafo out', 'poc']

    try:
        # the leading '@' used by the previous versions of the launcher is accepted and stripped
        CURTAILMENT['grid_limitation_point'] = ' '.join(str(CURTAILMENT['grid_limitation_point']).lower().replace('@', ' ').split())
    except:
        CURTAILMENT['grid_limitation_point'] = 'none'

    if CURTAILMENT['grid_limitation_point'] in ['', 'none']:
        CURTAILMENT['grid_limitation_point'] = None
    elif CURTAILMENT['grid_limitation_point'] not in valid_curtail_point:
        check_err_msg +='WARNING! Invalid grid limitation point. Valid points are: ' + str(valid_curtail_point) + ' or None. For this simulation it will be assumed None \n\n' 
        CURTAILMENT['grid_limitation_point'] = None

    # a curtailment applied downstream of parts that do not exist has no meaning
    if CURTAILMENT['grid_limitation_point'] in ['hv trafo out', 'poc'] and CONNECTION['hv_trafo_p_nom'] is None:
        check_err_msg +='WARNING! The grid limitation point was set downstream of the MV/HV transformer, which is not present. The curtailment has been moved to the entry of the HV transformer \n\n'
        CURTAILMENT['grid_limitation_point'] = 'hv trafo in'

    if CURTAILMENT['grid_limitation_point'] is None:
        CURTAILMENT['grid_curtail_power'] = None
    else:
        try:
            if not(CURTAILMENT['grid_curtail_power']>0):
                raise ValueError
        except:
            check_err_msg +='WARNING! Grid curtailment power - invalid or missing data. The simulation assumed the granted connection power \n\n'
            CURTAILMENT['grid_curtail_power'] = POWER['p_ac_connection']


    # check that Pdc, module, string lenght and inverter are compatible
    if not panerrorflag and not onderrorflag and not errorflag:
        Voc_string = MODULE['Voc'] * MMS['modules_per_string']*1.05  # a safety margin of 5% is included
        min_string_len = np.ceil((INVERTER['VMppMin'] / MODULE['Vmp'])*1.10)
        max_dc_ac_ratio = (MODULE['PNom'] * MMS['modules_per_string'] * INVERTER['NbInputs']) / INVERTER['PNomConv']/1000
        
        if Voc_string > INVERTER['VMPPMax']:
            max_string_len = np.floor(INVERTER['VMPPMax']/MODULE['Voc']/1.05)
            check_err_msg+='ERROR! Too many modules per string. Maximum number of module per string to match inverter limitation is: ' + str(max_string_len)+ '\n\n'
            errorflag = True
            
        if MMS['modules_per_string'] < min_string_len:
            check_err_msg+='ERROR! Too few modules per string. Minimum number of module per string to match inverter limitation is: ' + str(min_string_len)+ '\n\n'
            errorflag = True
            
        if MODULE['Isc'] > INVERTER['i_max_mppt']:
            check_err_msg+='ERROR! The module power is too high for ther inverter. Consider using a PV module with lower rated power or an invrter with higher MPPT current limit \n\n'
            errorflag = True
            
        try: 
            if (POWER['p_dc_peak']/POWER['p_ac_connection']) > max_dc_ac_ratio :
                check_err_msg+='ERROR! The desired DC/AC ratio is not theoretically achievable with the chosen module, string lenght and inverter. Reduce the DC/AC ratio or increase module per string\n\n'
                errorflag = True
        except:
            if max_dc_ac_ratio < 1.1 :
                check_err_msg+='WARNING! With the chosen module, string lenght and inverter the maximum DC/AC ratio achievable is 1.10 - Please review the design \n\n'

        # adjust Project peak to be equal to an integer number of strings
        POWER['p_dc_peak'] = round(POWER['p_dc_peak']*1000*1000/MODULE['PNom']/MMS['modules_per_string'])*MODULE['PNom']*MMS['modules_per_string']/1000/1000
        check_err_msg+='WARNING! The Project DC power was rounded to ' + str(POWER['p_dc_peak']) + ' MWp to match an integer number of strings \n\n'


    if errorflag:
        check_err_msg+='ERROR! The software encoutered some ERRORS and cannot continue the simulation. Please check your inputs \n'


    return POWER, LOCATION, METEO, MMS, GROUND, MODULE, INVERTER, TRAFO, LOSS, CONNECTION, CURTAILMENT, errorflag, check_err_msg



def pan_import_classic(file_path):
    """
    Reads a .PAN text file and extracts variables that compare in the file as `variable_name = value`.
    Infers the type of each variable (integer, float, boolean, string, list).
    Converts the value accordingly and returns a dictionary with variable names, their values, and types.
    
    Parameters
    ----------
    file_path : str or path object
        Name or path of a .pan file

    Returns
    -------
    MODULE : dict
        Dictionary containing PAN file variables, their values, and types.
    """

    MODULE = {}
    with open(file_path, 'r') as file:
        for line in file:
            # check if the ''End of PVObject pvIAMIAMProfile'' is in the line.
            # If so, stores the iam angle and iam loss in separated variables
            if 'End of TCubicProfile' in line:
                MODULE['iam_theta'] = np.zeros(shape = MODULE['NPtsEff'])
                MODULE['iam_loss'] = np.zeros(shape = MODULE['NPtsEff'])
                for n in np.arange(1, MODULE['NPtsEff']+1):
                    key = 'Point_' + str(n)
                    MODULE['iam_theta'][n-1] = MODULE.get(key)[0]
                    MODULE['iam_loss'][n-1] = MODULE.get(key)[1]
            
            elif '=' in line:  # Look for the presence of '='
                key, value = line.split('=', 1)  # Split into key and value
                value_cleaned = value.strip()  # Clean whitespace
                var_value = infer_and_convert_classic(value_cleaned)  # Get converted value
                MODULE[key.strip()] = var_value  # Store only the converted value

    # remove all unused variables from MODULE dict        
    keep_variable = ['Manufacturer', 'Model', 'PNom', 'Height', 'Width', 'NCelS', 'NCelP', 'NDiode', 'Isc', 'Voc', 'Imp', 'Vmp', 'muISC', 'Gamma', 'muGamma', 'RShunt', 'Rp_0', 'Rp_Exp', 'RSerie', 'Absorb', 'BifacialityFactor', 'iam_theta', 'iam_loss', 'SubModuleLayout', 'VRevDiode']
    MODULE = {key: value for key, value in MODULE.items() if key in keep_variable}

    return MODULE





def ond_import_classic(file_path):
    """
    Reads a .OND text file and extracts variables that compare in the file as `variable_name = value`.
    Infers the type of each variable (integer, float, boolean, string, list).
    Converts the value accordingly and returns a dictionary with variable names, their values, and types.
    
    Parameters
    ----------
    file_path : str or path object
        Name or path of a .pan file

    Returns
    -------
    INVERTER : dict
        Dictionary containing PAN file variables, their values, and types.
    """

    # first parse of the file
    # PROBLEM: Inverter efficiency variable have all same file name, so the parser will overwrite the first variable with the latter in the file having same name
    # SOLUTION: For each line, check if a section where inverter efficiency is start of finished
    with open(file_path, 'r') as file:
        INVERTER = {}
        V_eff_section = 0
        for line in file:
            
            # first check if the VNomEff is in the line.
            # If so, stores the voltage at which the efficiency is computed in separated variables
            if 'VNomEff' in line:
                key, VNomEff_list = line.split('=', 1)  # Split into key and value
                VNomEff_list = VNomEff_list.strip().split(',')
                for i, value in enumerate(VNomEff_list, start=1):
                    variable_name = 'v_eff'+ str(i)
                    if value.strip():      # exclude empty spaces from being converter to values
                        INVERTER[variable_name] = float(value)
            
            # then if the 'ProfilPIOV' is in the line, means that following lines will contain the efficiency values for a voltage
            # line with 'ProfilPIOV1'contain efficiency for lower voltage, 'ProfilPIOV2' for mid voltage and 'ProfilPIOV3' for very high voltage
            elif 'ProfilPIOV' in line:
                V_eff_section = V_eff_section + 1
            
            # if the 'End of TCubicProfile' is in the line and we were in a section where efficiency was outlined, it means the section is over and that voltage points can be resumed in a single variable
            # create the empt y array, sized according to the number of voltage provided
            elif 'End of TCubicProfile' in line and V_eff_section >0:
                INVERTER['p_in_v_eff' + str(V_eff_section)] = np.zeros(shape = INVERTER['NPtsEff'])
                INVERTER['p_out_v_eff' + str(V_eff_section)] = np.zeros(shape = INVERTER['NPtsEff'])
                
                for n in np.arange(1, INVERTER['NPtsEff']+1):
                    key = 'Point_' + str(n)
                    INVERTER['p_in_v_eff' + str(V_eff_section)][n-1] = INVERTER.get(key)[0]/1000   # converts W in kW
                    INVERTER['p_out_v_eff' + str(V_eff_section)][n-1] = INVERTER.get(key)[1]/1000  # converts W in kW
            
            # in all other case, just store the variable in the INVERTER dict
            elif '=' in line:  # Look for the presence of '='
                key, value = line.split('=', 1)  # Split into key and value
                value_cleaned = value.strip()  # Clean whitespace
                var_value = infer_and_convert_classic(value_cleaned)  # Get converted value
                INVERTER[key.strip()] = var_value  # Store only the converted value
    
    # calculate values from Inverter
    
    # remove all unused variables from INVERTER dict
    keep_variable = ['Manufacturer', 'Model', 'PNomConv', 'PMaxOUT','VMppMin', 'VMPPMax', 'IMaxDC', 'NbInputs', 'NbMPPT', 'VOutConv', 'TPMax','TPNom','TPLim1', 'TPLimAbs', 'PLim1', 'PLimAbs', 'v_eff1', 'v_eff2', 'v_eff3', 'p_in_v_eff1', 'p_in_v_eff2', 'p_in_v_eff3', 'p_out_v_eff1', 'p_out_v_eff2', 'p_out_v_eff3' ]
    INVERTER = {key: value for key, value in INVERTER.items() if key in keep_variable}

    # Return the variables dictionary
    return INVERTER





def infer_and_convert_classic(value: str):
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
    if "," in value:
        value_list = [infer_and_convert_classic(v.strip()) for v in value.split(",")]
        return value_list
    
    # 2. Boolean detection
    # Handles both True/False and Yes/No common conventions
    lower_value = value.lower()
    if lower_value in {"true", "false", "yes", "no"}:
        return lower_value in {"true", "yes"}
    
    # 3. Integer detection
    # Validates positive and negative whole numbers
    if value.isdigit() or (value.startswith("-") and value[1:].isdigit()):
        return int(value)
    
    # 4. Float detection
    # Tries to cast to float. Ensures it contains a decimal point or 
    # scientific notation ('e') to avoid overlapping with integers.
    try:
        float_value = float(value)
        if "." in value or "e" in lower_value:
            return float_value
    except ValueError:
        pass
    
    # 5. Default fallback
    # If all previous checks fail, return the original string
    return value





def post_processing_classic (RESULTS, energy_hourly, POWER, LOCATION, METEO, MMS, GROUND, MODULE, INVERTER, TRAFO, LOSS, CONNECTION, CURTAILMENT, errors):
    
    # 1 - IDENTIFY AVAILABLE FILE NAME
    current_date = datetime.now().strftime("%y%m%d")
    base_name = current_date + "_" + str("{:.2f}".format(POWER['p_dc_peak'])) + "MWp_" + str(POWER['p_ac_connection']) + "MVA_" + str(int(MODULE['PNom'])) + "Wp_" + str(int(INVERTER['PNomConv'])) + "kVA_" + MMS['mount_type'] + "_" + MMS['mod_configuration'] + "x" + str(MMS['modules_per_string']) +"_" + str(MMS['pitch']) + "m_pitch_" + str(MMS['tilt_angle']) + "°_tilt"
    extension = ".xlsx"
    
    results_dir = os.path.join(current_dir, 'Results')
    os.makedirs(results_dir, exist_ok=True)
    
    output_path = os.path.join(results_dir, base_name + extension)
    counter = 1
        
    # check if a file with same name exist. If so, it attach a _number to it            
    while os.path.exists(output_path):
        appendix = "_" + str(counter)
        new_name = base_name + appendix + extension
        output_path = os.path.join(results_dir, new_name)
        counter = counter + 1

    
    # 2 - CONVERT DICT INTO A DATAFRAME
    POWER = pd.DataFrame.from_dict(POWER, orient='index', columns=['Value'])
    LOCATION = pd.DataFrame.from_dict(LOCATION, orient='index', columns=['Value'])
    METEO = pd.DataFrame.from_dict(METEO, orient='index', columns=['Value'])
    MMS = pd.DataFrame.from_dict(MMS, orient='index', columns=['Value'])
    GROUND = pd.DataFrame.from_dict(GROUND, orient='index', columns=['Value'])
    MODULE = pd.DataFrame.from_dict(MODULE, orient='index', columns=['Value'])
    INVERTER = pd.DataFrame.from_dict(INVERTER, orient='index', columns=['Value'])
    TRAFO = pd.DataFrame.from_dict(TRAFO, orient='index', columns=['Value'])
    CONNECTION = pd.DataFrame.from_dict(CONNECTION, orient='index', columns=['Value'])
    CURTAILMENT = pd.DataFrame.from_dict(CURTAILMENT, orient='index', columns=['Value'])
    LOSS = pd.DataFrame.from_dict(LOSS, orient='index', columns=['Value'])
    
    # 3 - ADD UNIT COLUMN AND RENAME INDEX
    POWER_unit = {'p_dc_peak'       : 'MWp',
                  'p_ac_connection' : 'MWac'}
    POWER['Unit'] = POWER.index.map(POWER_unit)
    POWER = POWER.rename(index={'p_dc_peak'       : 'Plant peak power',
                                'p_ac_connection' : 'Connection AC power'})
    
    LOCATION_unit = {'latitude'  : 'deg',
                     'longitude' : 'deg'}
    LOCATION['Unit'] = LOCATION.index.map(LOCATION_unit)
    LOCATION = LOCATION.rename(index={'latitude' : 'Latitude',
                                      'longitude' : 'Longitude'})
    
    METEO = METEO.drop(index=['api_key'])   # remove API key for privacy
    METEO_unit = {'TMY_source'  : '-',
                  'TMY_file_name' : '-'}
    METEO['Unit'] = METEO.index.map(METEO_unit)
    METEO = METEO.rename(index={'TMY_source' : 'TMY Source',
                                'TMY file name (only for manual TMY import)' : '-'})
    
    MMS_unit = {'azis_azimuth'       : 'deg',
                'backtrack'          : '-',
                'ground_gap'         : 'm',
                'mod_configuration'  : '-',
                'modules_per_string' : '-',
                'mount_type'         : '-',
                'pitch'              : 'm',
                'tilt_angle'         : 'deg',
                'axis_azimuth'       : 'deg'}
    MMS['Unit'] = MMS.index.map(MMS_unit)
    MMS = MMS.rename(index={'azis_azimuth'       : 'Structure axis azimuth',
                            'backtrack'          : 'Backtracking active',
                            'ground_gap'         : 'Ground clearance @ max tilt',
                            'mod_configuration'  : 'Module layout configuration',
                            'modules_per_string' : 'Modules per string',
                            'mount_type'         : 'Structure type',
                            'pitch'              : 'Pitch',
                            'axis_azimuth'       : 'Axis azimuth',
                            'tilt_angle'         : 'Tilt angle'})
    
    GROUND_unit = {'albedo'        : '-',
                   'slope_azimuth' : 'deg',
                   'slope_tilt'    : 'deg'}
    GROUND['Unit'] = GROUND.index.map(GROUND_unit)
    GROUND = GROUND.rename(index={'albedo'        : 'Albedo',
                                  'slope_azimuth' : 'Slope azimuth',
                                  'slope_tilt'    : 'Slope tilt'})
    
    MODULE_unit = {'Manufacturer'  : '-',
                   'Model'         : '-',
                   'PAN_import'    : '-',
                   'PAN_file_name' : '-',
                   'PNom'          : 'Wp',
                   'Height'        : 'm',
                   'Width'         : 'm',
                   'NCelS'         : '-',
                   'NCelP'         : '-',
                   'NDiode'        : '-',
                   'Isc'           : 'A',
                   'Voc'           : 'V',
                   'Imp'           : 'A',
                   'Vmp'           : 'V',
                   'muISC'         : 'mAmps/Celsius',
                   'Gamma'         : '-',
                   'muGamma'       : '1/K',
                   'RShunt'        : 'ohm',
                   'Rp_0'          : 'ohm',
                   'Rp_Exp'        : '-',
                   'RSerie'        : 'ohm',
                   'Absorb'        : '-',
                   'BifacialityFactor': '-',
                   'iam_theta'     : 'deg',
                   'iam_loss'      : '-',
                   'VRevDiode'     : 'V',
                   'SubModuleLayout': '-'}
    MODULE['Unit'] = MODULE.index.map(MODULE_unit)
    MODULE = MODULE.rename(index={'Manufacturer'     : 'Manufacturer',
                                  'Model'            : 'Model',
                                  'PAN_import'       : 'Impoted user PAN file?',
                                  'PAN_file_name'    : 'PAN file name',
                                  'PNom'             : 'nominal power',
                                  'Height'           : 'Module eight',
                                  'Width'            : 'Module width',
                                  'NCelS'            : 'Number of cells in series',
                                  'NCelP'            : 'Number of cells in parallel',
                                  'NDiode'           : 'number of Bypass diodes',
                                  'Isc'              : 'Short circuit current',
                                  'Voc'              : 'Open circuit voltage',
                                  'Imp'              : 'max power operating current',
                                  'Vmp'              : 'max power operating voltage',
                                  'muISC'            : 'short-circuit current temperature coeff.',
                                  'Gamma'            : 'Diode ideality factor',
                                  'muGamma'          : 'Diode ideality factor temperature coeff.',
                                  'RShunt'           : 'Shunt resistance',
                                  'Rp_0'             : 'shunt resistance at zero irradiance',
                                  'Rp_Exp'           : 'Single diode exponential paramenter',
                                  'RSerie'           : 'series resistance at reference conditions',
                                  'Absorb'           : 'absorption coeff.',
                                  'BifacialityFactor': 'Bifaciality factor',
                                  'iam_theta'        : 'IAM loss reference angle',
                                  'iam_loss'         : 'IAM loss @ reference angle',
                                  'VRevDiode'        : 'bypass diode voltage drop',
                                  'SubModuleLayout'  : 'module internal configuration',})
    
    INVERTER_unit = {'Manufacturer'  : '-',
                     'Model'         : '-',
                     'OND_import'    : '-',
                     'OND_file_name' : '-',
                     'dc_ac_ratio'   : '-',
                     'PNomConv'      : 'kWac',
                     'PMaxOUT'       : 'kWac',
                     'VMppMin'       : 'V',
                     'VMPPMax'       : 'V',
                     'i_max_mppt'     : 'A',
                     'NbInputs'      : '-',
                     'NbMPPT'        : '-',
                     'VOutConv'      : 'V',
                     'inv_euro_eff'  : '%',
                     'IMaxDC'        : 'A',
                     'p_in_v_eff1'   : 'kWp',
                     'p_in_v_eff2'   : 'kWp',
                     'p_in_v_eff3'   : 'kWp',
                     'p_out_v_eff1'  : 'kWac',
                     'p_out_v_eff2'  : 'kWac',
                     'p_out_v_eff3'  : 'kWac',
                     'v_eff1'        : 'V',
                     'v_eff2'        : 'V',
                     'v_eff3'        : 'V'}
    INVERTER['Unit'] = INVERTER.index.map(INVERTER_unit)
    INVERTER = INVERTER.rename(index={'Manufacturer'  : 'Manufacturer',
                                      'Model'         : 'Model',
                                      'OND_import'    : 'Impoted user ONDfile?',
                                      'OND_file_name' : 'OND file name',
                                      'dc_ac_ratio'   : 'DC/AC ratio',
                                      'PNomConv'      : 'Nominal inverter AC power',
                                      'PMaxOUT'       : 'Maximum AC output power of inverter',
                                      'VMppMin'       : 'Minimum MPPT voltage',
                                      'VMPPMax'       : 'Maximum MPPT voltage',
                                      'i_max_mppt'     : 'Maximum current per MPPT',
                                      'NbInputs'      : 'Number of inverter inputs',
                                      'NbMPPT'        : 'Number of MPPT',
                                      'VOutConv'      : 'Inverter output voltage',
                                      'inv_euro_eff'  : 'Inverter EURO efficiency',
                                      'IMaxDC'        : 'Maximum total input current',
                                      'p_in_v_eff1'   : 'Efficiency curve 1: DC reference input power @ Veff_1',
                                      'p_in_v_eff2'   : 'Efficiency curve 2: DC reference input power @ Veff_2',
                                      'p_in_v_eff3'   : 'Efficiency curve 3: DC reference input power @ Veff_3',
                                      'p_out_v_eff1'  : 'Efficiency curve 1: AC reference input power @ Veff_1',
                                      'p_out_v_eff2'  : 'Efficiency curve 2: AC reference input power @ Veff_2',
                                      'p_out_v_eff3'  : 'Efficiency curve 3: AC reference input power @ Veff_3',
                                      'v_eff1'        : 'Efficiency curve 1: reference voltage Veff_1',
                                      'v_eff2'        : 'Efficiency curve 2: reference voltage Veff_2',
                                      'v_eff3'        : 'Efficiency curve 3: reference voltage Veff_3'})
    
    TRAFO_unit = {'p_trafo_max'      : 'MVA'}
    TRAFO['Unit'] = TRAFO.index.map(TRAFO_unit)
    TRAFO = TRAFO.rename(index={'p_trafo_max'      : 'LV/MV Trasformer max power'})

    LOSS_unit = {'soiling_loss'        : '%',
                 'rear_shade_factor'   : '%',
                 'thermal_loss_factor' : 'W/m2K',
                 'module_quality'      : '%',
                 'lid_loss'            : '%',
                 'first_year_degrad'   : '%',
                 'mismatch_loss_front' : '%',
                 'mismatch_loss_back'  : '%',
                 'dc_ohmic_loss_STC'   : '%',
                 'trafo_iron_loss'     : '%',
                 'trafo_copper_loss'   : '%',
                 'pv_availability'     : '%',
                 'grid_availability'   : '%',
                 'auxiliaries_loss'    : '%',
                 'ac_lv_cable_len'     : 'm',
                 'ac_lv_cable_sect'    : 'mm2',
                 'ac_mv_cable_len'     : 'm',
                 'ac_mv_cable_sect'    : 'mm2'}
    LOSS['Unit'] = LOSS.index.map(LOSS_unit)
    LOSS = LOSS.rename(index={'soiling_loss'        : 'Soiling loss',
                              'rear_shade_factor'   : 'Rear shade factor',
                              'thermal_loss_factor' : 'PV module thermal loss coeff.',
                              'module_quality'      : 'Module quality gain',
                              'lid_loss'            : 'LID loss',
                              'first_year_degrad'   : 'First year degradation',
                              'mismatch_loss_front' : 'Front mismatch loss',
                              'mismatch_loss_back'  : 'Back mismatch loss',
                              'dc_ohmic_loss_STC'   : 'DC ohmic wiring loss @STC',
                              'trafo_iron_loss'     : 'Iron loss trasformer',
                              'trafo_copper_loss'   : 'Copper loss Trasformer',
                              'pv_availability'     : 'Plant availability',
                              'grid_availability'   : 'Grid availability',
                              'auxiliaries_loss'    : 'Auxiliaries loss',
                              'ac_lv_cable_len'     : 'Average lenght of LV AC cables (from inverter to Trafo)',
                              'ac_lv_cable_sect'    : 'Section of LV AC cables',
                              'ac_mv_cable_len'     : 'Average lenght of MV AC cables (from trafo di collection cabin)',
                              'ac_mv_cable_sect'    : 'Section of MV AC cables'})

    CONNECTION_unit = {'mv_v_connection'         : 'kV',
                       'mv_connection_cable_len'  : 'km',
                       'mv_connection_cable_sect' : 'mm2',
                       'mv_connection_cable_nr'   : '-',
                       'hv_trafo_p_nom'           : 'MVA',
                       'hv_trafo_iron_loss'       : '%',
                       'hv_trafo_copper_loss'     : '%',
                       'hv_connection_cable_len'  : 'km'}
    CONNECTION['Unit'] = CONNECTION.index.map(CONNECTION_unit)
    CONNECTION = CONNECTION.rename(index={'mv_v_connection'         : 'MV connection voltage',
                                          'mv_connection_cable_len'  : 'Lenght of the MV connection cable, from collection cabin onwards',
                                          'mv_connection_cable_sect' : 'Section of the MV connection cable, from collection cabin onwards',
                                          'mv_connection_cable_nr'   : 'Number of MV connection cable, from collection cabin onwards',
                                          'hv_trafo_p_nom'           : 'MV/HV transformer rated power',
                                          'hv_trafo_iron_loss'       : 'Iron loss of the MV/HV transformer',
                                          'hv_trafo_copper_loss'     : 'Copper loss of the MV/HV transformer',
                                          'hv_connection_cable_len'  : 'Lenght of the HV connection line, from MV/HV transformer to PoC'})

    CURTAILMENT_unit = {'grid_limitation_point' : '-',
                        'grid_curtail_power'    : 'MVA'}
    CURTAILMENT['Unit'] = CURTAILMENT.index.map(CURTAILMENT_unit)
    CURTAILMENT = CURTAILMENT.rename(index={'grid_limitation_point' : 'Grid curtailment application point',
                                            'grid_curtail_power'    : 'Grid curtailment AC power'})

    error_list = [sentence.strip() for sentence in errors.split('\n\n') if sentence.strip()]
    
    # 4 - STORE VARIABLES IN AN EXCEL WORKBOOK
    OUTPUT = {'POWER'      : POWER,
              'LOCATION'   : LOCATION,
              'METEO'      : METEO,
              'MMS'        : MMS,
              'GROUND'     : GROUND,
              'MODULE'     : MODULE, 
              'INVERTER'   : INVERTER,
              'TRAFO'      : TRAFO,
              'CONNECTION' : CONNECTION,
              'CURTAILMENT': CURTAILMENT,
              'LOSS'       : LOSS,
              'RESULTS'    : RESULTS}
    
    workbook = Workbook()
    
    simulation_sheet = workbook.active
    simulation_sheet.title = 'Simulation'
    _build_simulation_sheet(simulation_sheet, OUTPUT)
    
    energy_sheet = workbook.create_sheet('Hourly energy')
    _build_energy_sheet(energy_sheet, energy_hourly)
    
    warnings_sheet = workbook.create_sheet('Warnings')
    _build_warnings_sheet(warnings_sheet, error_list)
    
    try:
        workbook.save(output_path)
    except PermissionError:
        raise PermissionError('Could not write "' + str(output_path) + '". The file is probably already open '
                              'in Excel, please close it and run the simulation again')
    
    print('\nResults exported to: ' + str(output_path))
    
    return output_path
            



# ==========================================================================
# ENGINE SUPPORT FUNCTIONS
# ==========================================================================
# The three functions below used to live at the bottom of PVyield_classic.py. They have been moved here so
# that the engine holds the simulation flow only, and every reusable piece of maths sits in this module
def get_irradiance_full_components (surface_tilt, surface_azimuth, solar_zenith, solar_azimuth, aoi, ghi, dhi, dni, dni_extra, airmass, albedo):
    # custom function to calculate full radiation components.
    # OUTPUT
    # components = dict with the following items:
        # 'poa_global'                  = global radiation on POA
        # 'poa_direct'                  = direct radiation on POA
        # 'poa_diffuse'                 = total diffuse radiation on POA (sky_diffuse + ground_diffuse)
        # 'poa_ground_diffuse'          = diffuse ground radiation on POA
        # 'poa_sky_diffuse'             = total sky diffuse radiation on POA (isotropic + circumsolar + horizon)
        # 'poa_sky_diffuse_isotropic'   = diffuse isotropic radiation (the one coming from sky dome)
        # 'poa_sky_diffuse_circumsolar' = diffuse  circumsolar radiation (the coming from area around the sun)
        # 'poa_sky_diffuse_horizon'     = diffuse horizon radiation (the one coming from 15° strip above horizon)
        
   
    
    #calculate full diffuse sky components (isotropic, circumsolar, horizon)
    diffuse_sky = irradiance.perez_driesse(surface_tilt = surface_tilt,
                                           surface_azimuth = surface_azimuth,
                                           solar_zenith = solar_zenith,
                                           solar_azimuth = solar_azimuth,
                                           dni = dni,
                                           dhi = dhi,
                                           dni_extra = dni_extra,
                                           airmass = airmass,
                                           return_components = True,
                                           )
    
    # set to 0 all negative horizon diffuse radiation contribution    
    diffuse_sky['poa_horizon'] = np.where(diffuse_sky['poa_horizon']<0, 0, diffuse_sky['poa_horizon'])
    
    #calculate ground diffuse component
    diffuse_ground = irradiance.get_ground_diffuse(surface_tilt = surface_tilt,
                                                   ghi = ghi, 
                                                   albedo = albedo, 
                                                   surface_type = None # If provided, overrightes albedo value. check Admissible value below
                                                   )
    # calculate direct component
    beam_component = np.maximum(dni * np.cos(np.radians(aoi)), 0)
    
    #allocate empty output variable to store all components
    components = pd.DataFrame()
    
    components['poa_global'] = beam_component + diffuse_ground + diffuse_sky['poa_sky_diffuse']
    components['poa_direct'] = beam_component
    components['poa_diffuse'] = diffuse_ground + diffuse_sky['poa_sky_diffuse']
    components['poa_ground_diffuse'] = diffuse_ground
    components['poa_sky_diffuse'] = diffuse_sky['poa_sky_diffuse']
    components['poa_sky_diffuse_isotropic'] = diffuse_sky['poa_isotropic']
    components['poa_sky_diffuse_circumsolar'] = diffuse_sky['poa_circumsolar']
    components['poa_sky_diffuse_horizon'] = diffuse_sky['poa_horizon']
    
    return components





def inv_Pout_from_Pin (v_eff1, v_eff2, v_eff3, p_in_v_eff1, p_in_v_eff2, p_in_v_eff3, p_out_v_eff1, p_out_v_eff2, p_out_v_eff3, p_dc_inv_in, v_dc_in):
    """
    Calculate the inverter AC output power given the input power, the MPPT input voltage and 3 efficiency curves data expressed as Pout = f(Vin, Pin)
    
    Parameters
    ----------
    v_eff1 : float, single value 
        inverter output voltage at which the efficiency is known for p_in_v_eff1 and p_out_v_eff1  [V]
    v_eff2 : float, single value 
        inverter output voltage at which the efficiency is known for p_in_v_eff2 and p_out_v_eff2  [V]
    v_eff3 : float, single value 
        inverter output voltage at which the efficiency is known for p_in_v_eff3 and p_out_v_eff3  [V]
    p_in_v_eff1 : numpy array
        Inverter input DC power at which the inverter efficiency is known for Voltage V_eff1 [Wp]
    p_in_v_eff2 : numpy array
        Inverter input DC power at which the inverter efficiency is known for Voltage V_eff2 [Wp]
    p_in_v_eff3 : numpy array
        Inverter input DC power at which the inverter efficiency is known for Voltage V_eff3 [Wp]
    p_out_v_eff1 : numpy array
        Inverter output AC power resulting from P_in_V_eff1 and V_eff1 [Wac]
    p_out_v_eff2 : numpy array
        Inverter output AC power resulting from P_in_V_eff2 and V_eff2 [Wac]
    p_out_v_eff3 : numpy array
        Inverter output AC power resulting from P_in_V_eff3 and V_eff3 [Wac]
    p_dc_inv_in : numpy array
        DC input power on inverter [Wp]
    v_dc_in : float, single value 
        inverter MPPT input voltage [V]

    Returns
    -------
    p_out_inv: numpy array
        AC output power on inverter [Wac]
    """  
    
    # Create interpolations functions for each voltage level P_out = f(P_in) for each voltage level
    # note that if the in the interpolation rules, we specify that if Pin is lower or higher that the boundaries expresses in the inverter efficiency curve, then P_out is clipped to minimum/max value
    interp_veff1 = interp1d(p_in_v_eff1, p_out_v_eff1, kind='quadratic', fill_value='extrapolate')
    interp_veff2 = interp1d(p_in_v_eff2, p_out_v_eff2, kind='quadratic', fill_value='extrapolate')
    interp_veff3 = interp1d(p_in_v_eff3, p_out_v_eff3, kind='quadratic', fill_value='extrapolate')
    
    # apply the interpolation functions above to the input power
    p_out1 = interp_veff1(p_dc_inv_in)
    p_out2 = interp_veff2(p_dc_inv_in)
    p_out3 = interp_veff3(p_dc_inv_in)
    
    # Stack the results for interpolation across voltages
    v_eff_values = [v_eff1, v_eff2, v_eff3] # create a single array with the 3 voltages
    p_out_values = np.vstack([p_out1, p_out2, p_out3]).T
    
    # Initialize an array for the final output power
    p_out_final = np.zeros(len(v_dc_in))
    
    # Interpolate for each target voltage
    for i, v_eff in enumerate(v_dc_in):
        if v_eff <= v_eff1:
            p_out_final[i] = p_out_values[i, 0]
        elif v_eff >= v_eff3:
            p_out_final[i] = p_out_values[i, 2]
        else:
            # Interpolate between voltages for the current index
            interp_voltages = interp1d(v_eff_values, p_out_values[i, :], kind='quadratic')
            p_out_final[i] = interp_voltages(v_eff)
    
    # Return the result as a pandas Series with the same index as the input Series
    return pd.Series(p_out_final, index=p_dc_inv_in.index)




def ac_ohm_loss (p_ac, v_p2p, leng, material, section):
    # calculate ac ohmic losses
    
    #INPUTS
    # p_ac = AC power [VA]
    # v_2p2 = phase-to-phase voltage [V]
    # lenght = cable lenght [m]
    # material = cable type. Either 'alu' or 'cu'
    # section = cable cross section [mm2]
    
    # OUTPUTS
    # p_loss = ohmic power loss between 0,1 [unitless]
    
    # calculate phase current
    i_p2g = p_ac/math.sqrt(3)/v_p2p
    
    # retuen resisivity as defined in PVsyst in ohm x m
    if material in ['cu', 'copper']:
        resistivity = 1.68E-8
    elif material in ['alu', 'aluminum']:
        resistivity = 2.82E-8
        
    resistance = resistivity / (section/np.power(100,3)) * leng
    
    ac_loss = 3* resistance * np.power(i_p2g,2) / p_ac      # cable loss for each simulation step    
    return ac_loss


def vf_row_ground_2d_weights (standard_angles, gcr, pvrow_height, pitch, npoints=100, n_pvrow=15, max_rows=4):
    """
    Computes how much each point of the ground contributes to the view of the ground from a row face.

    PVlib only exposes vf_row_ground_2d_integ, which is the view factor from the row to the WHOLE visible
    ground, already integrated. That is enough when the ground is uniformly lit, but not when it is not: the
    strip right under the rows is permanently shaded while the strip halfway between rows is bright, and a row
    does not look at the two equally. This function returns the missing piece, i.e. the SHAPE of that view, so
    that the brightness of the ground can be averaged with the proper weights instead of uniformly.

    The face is discretised into n_pvrow points. From each of them the contribution of a ground element is the
    2D differential view factor between a surface element and a strip

        dF = ( sin(theta_2) - sin(theta_1) ) / 2

    with the angles measured from the normal of the face. A ground element is discarded when the line of sight
    from the row point to it is intercepted by one of the neighbouring rows.

    The geometry follows the same convention as PVlib's vf_ground_sky_2d: rows infinitely long, regularly
    spaced by 'pitch', with the row rotated so that its RIGHT edge is the higher one. The result only depends
    on the geometry, so it is computed on a few structure tilt angles and interpolated afterwards.

    Front and rear are returned together because they share the whole geometry and differ only by the
    direction of the normal, so computing them in one pass costs almost nothing more than one of the two.

    NOTE: the weights are normalised to sum to one along the ground, i.e. they carry the SHAPE of the view and
    not its magnitude. The magnitude keeps coming from PVlib's vf_row_ground_2d_integ, so that any first order
    error in the absolute value of these weights cancels out.

    INPUT
    standard_angles : [deg] tilt angles of the STRUCTURE the geometry is computed on, from 0 (flat) to 90
    gcr             : [-] ground coverage ratio, row width divided by pitch
    pvrow_height    : [m] height of the centre of the PV row above the ground
    pitch           : [m] distance between the axis of two adjacent rows
    npoints         : [-] number of ground elements across one pitch
    n_pvrow         : [-] number of points the row face is discretised into
    max_rows        : [-] neighbouring rows considered on each side

    OUTPUT
    x_ground        : [-] centre of each ground element, as a fraction of the pitch, shape (npoints,)
    weights_front   : [-] contribution of each ground element to the view of the FRONT face,
                      shape (len(standard_angles), npoints), each row summing to one
    weights_back    : [-] same for the REAR face
    """

    angles = np.radians(np.atleast_1d(np.asarray(standard_angles, dtype=float)))
    half_width = gcr*pitch/2

    # the rows, as segments in the cross-axis plane, one per angle. The right edge is the higher one
    edge_right = np.stack([+half_width*np.cos(angles), pvrow_height + half_width*np.sin(angles)], axis=-1)
    edge_left = np.stack([-half_width*np.cos(angles), pvrow_height - half_width*np.sin(angles)], axis=-1)

    # points along the face, avoiding the very edges where the view factor is singular
    fraction = np.linspace(0.05, 0.95, n_pvrow)
    row_points = edge_left[:, np.newaxis, :] + fraction[np.newaxis, :, np.newaxis]*(edge_right - edge_left)[:, np.newaxis, :]

    # normal of each face. With the right edge higher, the front looks to the left and the rear to the right
    normal_front = np.stack([-np.sin(angles), np.cos(angles)], axis=-1)
    normals = np.stack([normal_front, -normal_front], axis=0)                  # (2, angles, 2), front then rear
    tangents = np.stack([-normals[..., 1], normals[..., 0]], axis=-1)

    # the ground, discretised across one pitch and centred below the row
    edges = np.linspace(-0.5, 0.5, npoints+1)*pitch
    x_centre = (edges[:-1] + edges[1:])/2

    # vectors from the row points to the two edges of each ground element.  shape (angles, n_pvrow, npoints)
    delta_x_1 = edges[np.newaxis, np.newaxis, :-1] - row_points[..., 0][..., np.newaxis]
    delta_x_2 = edges[np.newaxis, np.newaxis, 1:] - row_points[..., 0][..., np.newaxis]
    delta_y = -row_points[..., 1][..., np.newaxis]

    # the face index is prepended, so everything below carries the shape (2, angles, n_pvrow, npoints)
    normal_x = normals[..., 0][:, :, np.newaxis, np.newaxis]
    normal_y = normals[..., 1][:, :, np.newaxis, np.newaxis]
    tangent_x = tangents[..., 0][:, :, np.newaxis, np.newaxis]
    tangent_y = tangents[..., 1][:, :, np.newaxis, np.newaxis]

    angle_1 = np.arctan2(delta_x_1*tangent_x + delta_y*tangent_y, delta_x_1*normal_x + delta_y*normal_y)
    angle_2 = np.arctan2(delta_x_2*tangent_x + delta_y*tangent_y, delta_x_2*normal_x + delta_y*normal_y)

    # a ground element only counts when it lies in front of the face
    visible = (np.abs(angle_1) < np.pi/2) | (np.abs(angle_2) < np.pi/2)

    delta_vf = np.abs(np.sin(np.clip(angle_2, -np.pi/2, np.pi/2))
                      - np.sin(np.clip(angle_1, -np.pi/2, np.pi/2)))/2

    # lines of sight intercepted by the neighbouring rows, which are the same rows shifted by whole pitches
    shift = np.array([k*pitch for k in range(-max_rows, max_rows+1) if k != 0])
    blocker_1 = edge_left[:, np.newaxis, :] + np.stack([shift, np.zeros_like(shift)], axis=-1)[np.newaxis, :, :]
    blocker_2 = edge_right[:, np.newaxis, :] + np.stack([shift, np.zeros_like(shift)], axis=-1)[np.newaxis, :, :]

    p_x = row_points[..., 0][:, :, np.newaxis, np.newaxis]
    p_y = row_points[..., 1][:, :, np.newaxis, np.newaxis]
    c_x = x_centre[np.newaxis, np.newaxis, :, np.newaxis]
    c_y = 0.0
    q1_x = blocker_1[..., 0][:, np.newaxis, np.newaxis, :]
    q1_y = blocker_1[..., 1][:, np.newaxis, np.newaxis, :]
    q2_x = blocker_2[..., 0][:, np.newaxis, np.newaxis, :]
    q2_y = blocker_2[..., 1][:, np.newaxis, np.newaxis, :]

    def side (origin_x, origin_y, a_x, a_y, b_x, b_y):
        # cross product telling on which side of the segment origin-a the point b lies
        return (a_x-origin_x)*(b_y-origin_y) - (a_y-origin_y)*(b_x-origin_x)

    side_1 = side(p_x, p_y, c_x, c_y, q1_x, q1_y)
    side_2 = side(p_x, p_y, c_x, c_y, q2_x, q2_y)
    side_3 = side(q1_x, q1_y, q2_x, q2_y, p_x, p_y)
    side_4 = side(q1_x, q1_y, q2_x, q2_y, c_x, c_y)

    blocked = np.any((side_1*side_2 < 0) & (side_3*side_4 < 0), axis=-1)   # (angles, n_pvrow, npoints)

    # sum over the points of the face, then normalise along the ground
    weights = np.sum(np.where(visible & ~blocked[np.newaxis, ...], delta_vf, 0.0), axis=2)   # (2, angles, npoints)

    total = np.sum(weights, axis=-1, keepdims=True)
    weights = np.divide(weights, total, out=np.zeros_like(weights), where=total > 0)

    return x_centre/pitch, weights[0], weights[1]



# ==========================================================================
# EXCEL REPORT
# ==========================================================================
# Formatting constants, shared by the three sheet builders

_FONT_NAME = 'Arial'

_TITLE_FONT = Font(name=_FONT_NAME, size=14, bold=True, color='FFFFFF')
_TITLE_FILL = PatternFill('solid', fgColor='1F3864')

_SECTION_FONT = Font(name=_FONT_NAME, size=11, bold=True, color='FFFFFF')
_SECTION_FILL = PatternFill('solid', fgColor='2E75B6')

_HEADER_FONT = Font(name=_FONT_NAME, size=10, bold=True)
_HEADER_FILL = PatternFill('solid', fgColor='D9E2F3')

_BODY_FONT = Font(name=_FONT_NAME, size=10)
_BOLD_FONT = Font(name=_FONT_NAME, size=10, bold=True)
_NOTE_FONT = Font(name=_FONT_NAME, size=9, italic=True, color='595959')

_HIGHLIGHT_FILL = PatternFill('solid', fgColor='FFF2CC')
_THIN_BORDER = Border(bottom=Side(style='thin', color='BFBFBF'))

_NUM_FORMAT = '#,##0.00'
_ENERGY_FORMAT = '#,##0.0'

# rows of the results that close a block, highlighted in the report
_HIGHLIGHT_ROWS = ['PR PV only', 'Yield PV only',
                   'PR PV + MV connection + PV avail.', 'Yield PV + MV connection + PV avail.',
                   'PR PV + MV&HV connection + grid avail.', 'Yield PV + MV&HV connection + grid avail.']


def _write_title (sheet, row, text, n_columns):
    # writes a full width title bar and returns the next free row
    sheet.cell(row=row, column=1, value=text)
    for col in range(1, n_columns+1):
        sheet.cell(row=row, column=col).fill = _TITLE_FILL
        sheet.cell(row=row, column=col).font = _TITLE_FONT
    sheet.row_dimensions[row].height = 22
    return row + 1


def _write_section (sheet, row, text, n_columns):
    # writes a section bar and returns the next free row
    sheet.cell(row=row, column=1, value=text)
    for col in range(1, n_columns+1):
        sheet.cell(row=row, column=col).fill = _SECTION_FILL
        sheet.cell(row=row, column=col).font = _SECTION_FONT
    return row + 1


def _write_header (sheet, row, labels):
    # writes a table header and returns the next free row
    for col, label in enumerate(labels, start=1):
        cell = sheet.cell(row=row, column=col, value=label)
        cell.font = _HEADER_FONT
        cell.fill = _HEADER_FILL
        cell.border = _THIN_BORDER
        cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
    return row + 1


def _format_value (value):
    # converts a value read from the input dictionaries into something Excel can store
    if value is None:
        return 'not used'
    if isinstance(value, bool):
        return 'Yes' if value else 'No'
    if isinstance(value, (tuple, list, np.ndarray)):
        values = list(np.ravel(value))
        if len(values) > 12:
            return str(len(values)) + ' points (see source file)'
        return ', '.join(('%g' % v) if isinstance(v, (int, float, np.number)) else str(v) for v in values)
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    return value


def _build_simulation_sheet (sheet, OUTPUT):
    """
    Fills the first sheet with the simulation inputs and the loss chain.

    INPUT
    sheet  : destination openpyxl worksheet
    OUTPUT : dictionary of DataFrames, one per input group plus 'RESULTS', each holding a 'Value' and a
             'Unit' column, exactly as built by post_processing_classic
    """

    n_columns = 4

    row = _write_title(sheet, 1, 'PVyield classic  -  simulation inputs and results', n_columns)
    sheet.cell(row=row, column=1, value='Generated on ' + datetime.now().strftime('%Y-%m-%d %H:%M')).font = _NOTE_FONT
    row += 2

    # ---- inputs ----
    row = _write_section(sheet, row, 'SIMULATION INPUTS', n_columns)
    row = _write_header(sheet, row, ['Group', 'Parameter', 'Unit', 'Value'])

    for group_name, table in OUTPUT.items():
        if group_name == 'RESULTS':
            continue

        for index, table_row in table.iterrows():
            sheet.cell(row=row, column=1, value=group_name).font = _NOTE_FONT
            sheet.cell(row=row, column=2, value=str(index)).font = _BODY_FONT

            cell = sheet.cell(row=row, column=3, value=table_row['Unit'])
            cell.font = _BODY_FONT
            cell.alignment = Alignment(horizontal='center')

            cell = sheet.cell(row=row, column=4, value=_format_value(table_row['Value']))
            cell.font = _BODY_FONT
            if isinstance(cell.value, float):
                cell.number_format = _NUM_FORMAT

            row += 1

        row += 1   # a blank line separates the groups

    row += 1

    # ---- results ----
    row = _write_section(sheet, row, 'SIMULATION RESULTS  -  loss chain', n_columns)
    row = _write_header(sheet, row, ['Loss / gain', 'Unit', 'Value'])

    results = OUTPUT['RESULTS']

    for index, table_row in results.iterrows():

        # the spacers of the loss chain become blank separator lines
        if str(index).strip('- ') == '':
            row += 1
            continue

        is_highlight = index in _HIGHLIGHT_ROWS
        font = _BOLD_FONT if is_highlight else _BODY_FONT

        sheet.cell(row=row, column=1, value=str(index)).font = font

        cell = sheet.cell(row=row, column=2, value=table_row['Unit'])
        cell.font = font
        cell.alignment = Alignment(horizontal='center')

        value = table_row['Value']
        cell = sheet.cell(row=row, column=3, value=None if pd.isna(value) else float(value))
        cell.font = font
        cell.number_format = _NUM_FORMAT
        cell.alignment = Alignment(horizontal='center')

        if is_highlight:
            for col in range(1, 4):
                sheet.cell(row=row, column=col).fill = _HIGHLIGHT_FILL

        row += 1

    row += 1
    sheet.cell(row=row, column=1, value='Note: a negative value is a loss, a positive value is a gain. The '
                                        'Performance Ratio is the ratio between the Net Yield and the global '
                                        'irradiation on the collector plane, the same definition PVsyst '
                                        'uses, and not the product of the rows above.').font = _NOTE_FONT

    sheet.column_dimensions['A'].width = 46
    sheet.column_dimensions['B'].width = 58
    sheet.column_dimensions['C'].width = 12
    sheet.column_dimensions['D'].width = 22

    sheet.freeze_panes = 'A4'


def _build_energy_sheet (sheet, energy_hourly):
    """
    Fills the second sheet with the hourly net energy and a monthly summary built with live formulas.

    INPUT
    sheet         : destination openpyxl worksheet
    energy_hourly : [kWh] hourly net energy at the Point of Connection, as a pandas Series
    """

    _write_title(sheet, 1, 'PVyield classic  -  hourly net energy at the Point of Connection', 6)

    header_row = 3
    _write_header(sheet, header_row, ['Timestamp (local time)', 'Month', 'Net energy [kWh]'])

    # Excel cannot store a timezone aware timestamp, so the timezone is removed. The values are already
    # expressed in the local time of the site
    timestamps = energy_hourly.index
    if timestamps.tz is not None:
        timestamps = timestamps.tz_localize(None)

    first_data_row = header_row + 1

    for index, (timestamp, value) in enumerate(zip(timestamps, energy_hourly.to_numpy())):
        current_row = first_data_row + index

        cell = sheet.cell(row=current_row, column=1, value=timestamp.to_pydatetime())
        cell.font = _BODY_FONT
        cell.number_format = 'DD/MM/YYYY HH:MM'

        cell = sheet.cell(row=current_row, column=2, value=int(timestamp.month))
        cell.font = _BODY_FONT
        cell.alignment = Alignment(horizontal='center')

        cell = sheet.cell(row=current_row, column=3, value=float(np.nan_to_num(value)))
        cell.font = _BODY_FONT
        cell.number_format = _ENERGY_FORMAT

    last_data_row = first_data_row + len(energy_hourly) - 1

    # monthly summary, built with live formulas so that it follows any edit made to the hourly values
    month_names = ['January', 'February', 'March', 'April', 'May', 'June',
                   'July', 'August', 'September', 'October', 'November', 'December']

    summary_column = 5   # column E

    for offset, label in enumerate(['Month', 'Net energy [MWh]']):
        cell = sheet.cell(row=header_row, column=summary_column+offset, value=label)
        cell.font = _HEADER_FONT
        cell.fill = _HEADER_FILL
        cell.border = _THIN_BORDER
        cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)

    for month_index, month_name in enumerate(month_names):
        current_row = first_data_row + month_index

        cell = sheet.cell(row=current_row, column=summary_column, value=month_name)
        cell.font = _BODY_FONT

        formula = ('=SUMIFS($C$' + str(first_data_row) + ':$C$' + str(last_data_row) + ','
                   '$B$' + str(first_data_row) + ':$B$' + str(last_data_row) + ',' + str(month_index+1) + ')/1000')
        cell = sheet.cell(row=current_row, column=summary_column+1, value=formula)
        cell.font = _BODY_FONT
        cell.number_format = _ENERGY_FORMAT

    total_row = first_data_row + 12
    cell = sheet.cell(row=total_row, column=summary_column, value='TOTAL')
    cell.font = _BOLD_FONT
    cell.fill = _HIGHLIGHT_FILL

    total_letter = get_column_letter(summary_column+1)
    cell = sheet.cell(row=total_row, column=summary_column+1,
                      value='=SUM(' + total_letter + str(first_data_row) + ':' + total_letter + str(total_row-1) + ')')
    cell.font = _BOLD_FONT
    cell.fill = _HIGHLIGHT_FILL
    cell.number_format = _ENERGY_FORMAT

    sheet.column_dimensions['A'].width = 22
    sheet.column_dimensions['B'].width = 8
    sheet.column_dimensions['C'].width = 16
    sheet.column_dimensions['D'].width = 4
    sheet.column_dimensions['E'].width = 14
    sheet.column_dimensions['F'].width = 18

    sheet.freeze_panes = 'A' + str(first_data_row)


def _build_warnings_sheet (sheet, error_list):
    """
    Fills the third sheet with the errors and warnings raised while the simulation was running, one per row.

    INPUT
    sheet      : destination openpyxl worksheet
    error_list : list of message strings, already split and stripped
    """

    _write_title(sheet, 1, 'PVyield classic  -  errors and warnings raised during the simulation', 3)

    if not error_list:
        cell = sheet.cell(row=3, column=1, value='No warning was raised: every input was read as provided, '
                                                 'with no standard value assumed in its place.')
        cell.font = _BODY_FONT
        sheet.column_dimensions['A'].width = 100
        return

    sheet.cell(row=2, column=1, value=str(len(error_list)) + ' message(s) raised. A WARNING reports an input '
                                      'that was missing or inconsistent and has been replaced by a standard '
                                      'value, or a result worth checking. An ERROR stops the '
                                      'simulation.').font = _NOTE_FONT

    header_row = 4
    _write_header(sheet, header_row, ['#', 'Type', 'Message'])

    for index, message in enumerate(error_list):
        current_row = header_row + 1 + index

        cell = sheet.cell(row=current_row, column=1, value=index+1)
        cell.font = _BODY_FONT
        cell.alignment = Alignment(horizontal='center', vertical='top')

        message_type = 'ERROR' if message.upper().startswith('ERROR') else 'WARNING'
        cell = sheet.cell(row=current_row, column=2, value=message_type)
        cell.font = _BODY_FONT
        cell.alignment = Alignment(horizontal='center', vertical='top')

        cell = sheet.cell(row=current_row, column=3, value=message)
        cell.font = _BODY_FONT
        cell.alignment = Alignment(vertical='top', wrap_text=True)

    sheet.column_dimensions['A'].width = 6
    sheet.column_dimensions['B'].width = 12
    sheet.column_dimensions['C'].width = 120

    sheet.freeze_panes = 'A' + str(header_row+1)
