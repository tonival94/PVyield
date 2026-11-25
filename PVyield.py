# -*- coding: utf-8 -*-
"""
PVyield v15

Multi-field PV energy yield simulation engine
@author: Antonio Vallone (-T O N I-)

Libraries requirements:
> PVlib == 0.16.1
> pandas == 3.0.5
> scipy == 1.18.0
> timezonefinder == 8.5.2
> tzdata == 2026c

# improvements for next versions:
# 1 - Add various source for horizon import (e.g. solargis, solcast, etc...) - Waiting for PVlib to release functions
# 2 - Add various API functionality to automatically import TMY from their server (missing: SolarGIS, SolarAnywhere)  - Waiting for PVlib to release functions
# 3 - Add functions to read third party tmy (e.g. solcast). Need a TMY to check how data is stored
# 4 - Add Snow losses calculation (PVlib already includes functions to do that) - snow data can be imported from Solcast API
# 5 - include low irradiation correction in Single diode model. This will fix irradiation losses
# 6 - Calculate partial shading loss accurately by performing electric calculation of the series/parallel substrings in the module
# 7 - Add function to automatically calculate MV connection cable section and number (MV_Connection.auto_sizing)
"""

import csv
import warnings
import numpy as np
import pandas as pd
import math
from pathlib import Path
from datetime import datetime
from scipy.interpolate import interp1d
from timezonefinder import TimezoneFinder
import pytz

from pvlib import location
from pvlib import atmosphere
from pvlib import irradiance
from pvlib import iotools
from pvlib import bifacial
from pvlib.bifacial import utils
from pvlib import pvsystem
from pvlib import iam
from pvlib import shading
from pvlib import temperature
from pvlib import solarposition
from pvlib import tools
from pvlib import tracking
from pvlib import transformer

from PVyield_classes import Constants, Results

_CURRENT_DIR = Path(__file__).parent.resolve()   # get current file directory


def PVyield(pv_project):
    """
    Runs the PV energy yield simulation of a whole PV project.

    Args:
        pv_project (PV_Project):   PV_Project class object, already validated through input_check()

    Returns:
        RESULTS (pd.DataFrame):    loss chain and Performance Ratio for each PV field and for the whole project
        energy_hourly (pd.Series): [kWh] net energy delivered by the project at the Point of Connection,
                                   for each of the 8760 hours of the year
    """

    # =================================================================================================
    # 1 - PREALLOCATE RESULTS
    # =================================================================================================

    RESULTS = Results.assign_results(pv_project)

    # Bookkeeping used to build the 'PV project' column of the loss chain.
    # For every loss row we accumulate, across all the fields, the energy entering and the energy leaving that stage.
    # The project loss is then the ratio between the two sums, which is exactly an energy weighted average of the
    # individual field losses and keeps the product of the chain consistent with the project Net Yield.
    energy_in = {}
    energy_out = {}

    def record_loss(row, field_name, e_in, e_out):
        """
        Stores one loss row for one field and accumulates the energies needed for the project column.

        Args:
            row (str):         name of the RESULTS row
            field_name (str):  name of the PV field
            e_in (float):      energy entering the stage, in Wh equivalent for the whole field
            e_out (float):     energy leaving the stage, in Wh equivalent for the whole field
        """
        RESULTS.loc[row, field_name] = e_out/e_in - 1 if e_in else 0
        energy_in[row] = energy_in.get(row, 0) + e_in
        energy_out[row] = energy_out.get(row, 0) + e_out


    # =================================================================================================
    # 2 - GENERATE TIMESTEPS FOR SIMULATION
    # =================================================================================================

    latitude = pv_project.coordinates.latitude
    longitude = pv_project.coordinates.longitude

    # identify the timezone from the project coordinates
    timezone_name = TimezoneFinder().timezone_at(lat = latitude,
                                                 lng = longitude,
                                                 )

    timezone = pytz.timezone(timezone_name)

    # get the UTC offset, needed to realign the TMY downloaded from online databases (which is delivered in UTC)
    UTC_offset = timezone.utcoffset(datetime.now()).total_seconds() / 3600   # offset in hours

    local_time = pd.date_range(start = f'{Constants.sim_year}-01-01 00:00',
                               end = f'{Constants.sim_year}-12-31 23:00',
                               freq = Constants.freq,
                               tz = timezone,
                               )


    # =================================================================================================
    # 3 - SET LOCATION
    # =================================================================================================

    try:
        altitude = location.lookup_altitude(latitude = latitude,
                                            longitude = longitude,
                                            )
    except Exception:
        warnings.warn(f'Not able to retrieve altitude for the given location. '
                      f'The simulation assumed an altitude of {Constants.default_altitude} m above sea level')
        altitude = Constants.default_altitude

    site_location = location.Location(latitude = latitude,
                                      longitude = longitude,
                                      tz = timezone,
                                      name = None,
                                      altitude = altitude,
                                      )


    # =================================================================================================
    # 4 - IMPORT RADIATION DATA
    # =================================================================================================

    TMY_source = pv_project.tmy.TMY_source
    api_key = pv_project.tmy.api_key
    TMY_file_name = pv_project.tmy.TMY_file_name

    # 4.1 - IMPORT TMY DATA FROM THE SELECTED SOURCE
    try:
        if TMY_source == 'pvgis':
            imported_tmy, *_ = iotools.get_pvgis_tmy(latitude = latitude,
                                                     longitude = longitude,
                                                     outputformat = 'csv',
                                                     usehorizon = False,
                                                     userhorizon = None,
                                                     startyear = 2005,
                                                     endyear = 2015,
                                                     map_variables = True,
                                                     url = 'https://re.jrc.ec.europa.eu/api/',
                                                     timeout = 30,
                                                     )

        elif TMY_source == 'solcast api':
            # To test the API with no paid subscription, register to Solcast and then use coordinates:
            # latitude =  -33.856784
            # longitude = 151.215297
            imported_tmy, *_ = iotools.get_solcast_tmy(latitude = latitude,
                                                       longitude = longitude,
                                                       api_key = api_key,
                                                       map_variables = True,
                                                       output_parameters = ['ghi','dni','dhi','air_temp','surface_pressure','relative_humidity','wind_speed_10m','albedo'],
                                                       # Other output can be specified as specified in https://docs.solcast.com.au/#output-parameter-definitions
                                                       )

        elif TMY_source == 'meteonorm api':
            # To test the API with no paid subscription, use:
            # latitude =  50
            # longitude = 10
            # api_key = 'demo0000-0000-0000-0000-000000000000'
            imported_tmy, *_ = iotools.get_meteonorm_tmy(latitude = latitude,
                                                         longitude = longitude,
                                                         api_key = api_key,
                                                         parameters = 'all',
                                                         time_step = '1h',
                                                         horizon = 'flat',
                                                         map_variables = True,
                                                         # Other options are described in https://docs.meteonorm.com/api/climate#climate-tmy
                                                         )

            # the PVlib map_variables of Meteonorm does not cover the air pressure, so it keeps its original Meteonorm name and its original unit hPa.
            # It is renamed and converted to Pa here. NOTE: The parameter name is looked up rather than hardcoded, because it has changed across
            # Meteonorm API versions. If no pressure column is found, the standard value is used further below
            for meteonorm_pressure_name in ['surface_air_pressure', 'air_pressure', 'pressure', 'atmospheric_pressure']:
                if meteonorm_pressure_name in imported_tmy:
                    imported_tmy = imported_tmy.rename(columns={meteonorm_pressure_name: 'pressure'})
                    imported_tmy['pressure'] = imported_tmy['pressure'] * 100
                    break

        if TMY_source in ['pvgis', 'solcast api', 'meteonorm api']:

            # the TMY must hold exactly one full non leap year at the simulation timestep, otherwise the index
            # assignment below would fail with an opaque length mismatch error
            if len(imported_tmy) != len(local_time):
                raise RuntimeError(f'The TMY imported from "{TMY_source}" holds {len(imported_tmy)} records, while '
                                   f'{len(local_time)} were expected (one full non leap year at a {Constants.freq} timestep). '
                                   f'Please check the TMY source settings')

            # The online databases above deliver the TMY in UTC. The simulation runs in local time,
            # therefore the data are rolled by the UTC offset of the site.
            # NOTE: the roll is applied column by column so that the dtype of each column is preserved
            imported_tmy = imported_tmy.apply(lambda column: np.roll(column.to_numpy(), int(round(UTC_offset))))
            imported_tmy.index = local_time

    except Exception as err:
        raise RuntimeError(f'Not able to get radiation data from the online database "{TMY_source}". '
                           f'Please check your internet connection, your API key and the project coordinates. Original error: {err}')

    if TMY_source == 'solargis tmy':
        TMY_path = _CURRENT_DIR / 'TMY' / TMY_file_name

        try:
            # find the row where the actual TMY data starts
            with open(TMY_path) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=';')
                line_count = 1
                for row in csv_reader:
                    if row[0] == '#Data:':
                        break
                    else:
                        line_count = line_count+1

            # read the TMY file
            imported_tmy = pd.read_csv(TMY_path,
                                       sep = ';',
                                       header = line_count,
                                       )
        except Exception as err:
            raise RuntimeError(f'Not able to read the SolarGIS TMY file provided. Please make sure the file "{TMY_file_name}" '
                               f'exists inside the "TMY" folder and that it is a legit SolarGIS export. Original error: {err}')

        # SolarGIS TMY files are already in local time, no need to shift the data to match our timesteps
        imported_tmy.index = local_time

        # remap the variables to match the PVlib nomenclature
        imported_tmy = imported_tmy.rename(columns={'GHI':  'ghi',
                                                    'DIF':  'dhi',
                                                    'DNI':  'dni',
                                                    'TEMP': 'temp_air',
                                                    'AP':   'pressure',
                                                    'RH':   'relative_humidity',
                                                    'WS':   'wind_speed'}
                                           )

        # SolarGIS TMY reports the air pressure (AP) in hPa, while every PVlib function expects it in Pa.
        if 'pressure' in imported_tmy:
            imported_tmy['pressure'] = imported_tmy['pressure'] * 100

    # 4.2 - UNPACK DATA INTO LOCAL VARIABLES
    ghi = imported_tmy['ghi']
    dhi = imported_tmy['dhi']
    dni = imported_tmy['dni']
    temp_air = imported_tmy['temp_air']
    wind_speed = imported_tmy['wind_speed']

    # check whether the optional variables have been imported (not all the TMY sources carry them)
    if 'albedo' in imported_tmy:
        # when the TMY carries an albedo timeseries, it overrides the albedo declared in the Ground class of every field
        tmy_albedo = imported_tmy['albedo']
        warnings.warn(f'The TMY imported from "{TMY_source}" carries albedo data. '
                      f'The albedo declared in the Ground class of every PV field has been overridden by the TMY values')
    else:
        tmy_albedo = None

    if 'pressure' in imported_tmy:
        pressure = imported_tmy['pressure']

        # Sanity check on the unit. Every source is converted to Pa where needed, so a pressure of a few hundred
        # hPa here means a provider changed its unit, or a new source was added without its own conversion.
        # Sea level pressure is about 101325 Pa and it stays well above 50000 Pa even on the highest plateaus,
        # so anything below that threshold can only be a unit mismatch
        if np.nanmedian(pressure) < 50000:
            warnings.warn(f'The air pressure imported from "{TMY_source}" has a median of {np.nanmedian(pressure):.0f}, '
                          f'which is far below the {Constants.default_pressure} Pa expected at sea level. The data are '
                          f'probably expressed in hPa rather than Pa. The simulation used the standard pressure instead')
            pressure = pd.Series(Constants.default_pressure, index=local_time, name='air pressure')
    else:
        # if the TMY carries no pressure, the average air pressure is used instead
        pressure = pd.Series(Constants.default_pressure, index=local_time, name='air pressure')

    if 'relative_humidity' in imported_tmy:
        relative_humidity = imported_tmy['relative_humidity']   # at the moment not used

    # GHI is identical for every field of the project, being the project a single location.
    # Only the value columns are written, the 'Units' column holds strings and must not be touched
    RESULTS.loc[Results.GHI, Results.value_columns(RESULTS)] = np.nansum(ghi)/1000


    # =================================================================================================
    # 5 - COMPUTE TIME SHIFT
    # =================================================================================================
    """
    In the PVlib convention, the time label always refers to the beginning of the record, and concerns the accumulation up to the next record
    If times = 7:00, means the radiation data are the sum of radiation from 7:00 to 8:00, and solar positions are computed at 7:00
    Now, with imported data files, there is no general rule.
    Each record usually holds a time label, which possibly doesn't correspond to the PVlib definition.
    Moreover, this time label may represent data measured (averaged) over a time period which doesn't correspond to the true Legal Time "round" hours
    We need to align the imported TMY with our timesteps. This is done by trying to shift the imported data and see what shift gives a best match with clear sky irradiation
    """

    # for TMY source SolarGIS TMY and Solcast API the timeshift is known
    if TMY_source in ['SolarGIS TMY', 'Solcast API']:
        time_shift = pd.Timedelta(minutes = 30)

    else:
        # create a vector holding every timeshift that will be tried
        scan_start, scan_stop = Constants.time_shift_scan
        err = pd.Series(index = np.arange(scan_start, scan_stop+1, 1), dtype=float)

        # cycle through the timeshift range and compute the Mean Square Error for each tentative shift
        for time_shift in np.arange(scan_start, scan_stop, 1):
            time_delta = pd.Timedelta(minutes=int(time_shift))

            clear_sky = site_location.get_clearsky(times = local_time + time_delta,
                                                   model = 'ineichen',
                                                   )
            # assign to clear_sky the same index as ghi to allow for numpy calculation with ghi
            # NOTE: clear_sky is computed with the timeshift applied, we simply reassign the original index so that we can perform operations with the ghi variable
            clear_sky.index = local_time

            # compute the Mean Square Error and store it
            MSE = np.power(clear_sky['ghi']-imported_tmy['ghi'], 2).mean()
            err.loc[time_shift] = MSE

        time_shift = err.idxmin()                        # find the time_shift associated to the minimum MSE
        time_shift = pd.Timedelta(minutes=int(time_shift))    # convert the timeshift in a pandas timedelta

    # compute the new shifted timesteps that will be used for the solar geometry calculations
    local_time_calc = local_time + time_shift


    # =================================================================================================
    # 6 - CALCULATE SOLAR GEOMETRY
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

    solar_position.index = local_time   # change to local time for better readability

    # compute the extraterrestrial normal irradiance
    dni_extra = irradiance.get_extra_radiation(datetime_or_doy = local_time_calc,
                                               solar_constant = Constants.solar_constant,
                                               method = 'spencer',
                                               epoch_year = datetime.now().year,
                                               )

    dni_extra.index = local_time   # change to local time for better readability

    # compute the relative airmass
    airmass_relative = atmosphere.get_relative_airmass(zenith = solar_position['apparent_zenith'],
                                                       model = 'gueymard2003',
                                                       )


    # =================================================================================================
    # 7 - HORIZON DATA & LOSSES ON THE HORIZONTAL COMPONENTS
    # =================================================================================================
    # The horizon profile depends on the location only, therefore the horizon losses on the horizontal
    # irradiance components (DNI and DHI) are computed once for the whole project, outside the field loop.
    # The horizon loss ON THE PLANE OF ARRAY, on the contrary, depends on the orientation of the structure
    # and is therefore computed field by field.

    # 7.1 - GET HORIZON DATA
    # get the horizon profile from PVGIS. The result is a tuple
    try:
        hor_data = iotools.get_pvgis_horizon(latitude = latitude,
                                             longitude = longitude,
                                             url = 'https://re.jrc.ec.europa.eu/api/',
                                             )

    # if the horizon cannot be retrieved, a flat horizon is assumed
    except Exception:
        hor_data = [pd.Series(data = [0,0], index = [0,180], name = 'horizon_elevation'), 'Not able to retrieve horizon data from PVGIS. The simulation assumed flat horizon']
        warnings.warn('Could not get horizon elevation data from PVGIS. The simulation assumed a flat horizon')

    # extract the data from the tuple and create a DataFrame to store the future horizon data
    hor_data = pd.DataFrame(hor_data[0], index = hor_data[0].index)

    # 7.2 - GET IDEAL SUNRISE AND SUNSET TIMES
    # for each day, compute the ideal sunrise and sunset time, corresponding to a flat horizon
    day_times = pd.date_range(start = f'{Constants.sim_year}-01-01 00:00',
                              end = f'{Constants.sim_year}-12-31 23:59',
                              freq = '1D',
                              tz = timezone,
                              )

    sun_rise_set_year = solarposition.sun_rise_set_transit_spa(times = day_times + pd.Timedelta(days=0.5),   # by specifying the hour as noon (+0.5 days respect to the original timestep), the function returns the correct day for sunrise and sunset. If midnight is passed, the sunrise is calculated on the previous day
                                                               longitude = longitude,
                                                               latitude = latitude,
                                                               )

    # rename the DataFrame columns for better readability and change the index to the one matching the rest of the PVlib calculations
    sun_rise_set_year = sun_rise_set_year.rename(columns={'sunrise':'ideal_sunrise', 'sunset':'ideal_sunset'})
    sun_rise_set_year = sun_rise_set_year.drop('transit', axis = 1)
    sun_rise_set_year.index = day_times

    # 7.3 - GET REAL SUNRISE AND SUNSET TIMES
    # the real sunrise and sunset are computed from the imported horizon data
    for day in day_times:

        precise_time = pd.date_range(start = day,
                                     end = day + pd.Timedelta(hours=23, minutes=59),
                                     freq = '1min',
                                     tz = timezone,
                                     )

        # get the current hour of precise_time. To avoid issues with DST (Daylight Saving Time) we must first delocalise the data, then set minutes and seconds to zero, then relocalise the data
        current_hour = precise_time.tz_localize(None).floor('h').tz_localize(timezone, ambiguous = False)

        # compute the precise sun position for every minute. This will be used to get the precise time when the sun crosses the horizon
        sun_positions_rise_set = solarposition.get_solarposition(time = precise_time,
                                                                 latitude = latitude,
                                                                 longitude = longitude,
                                                                 altitude = altitude,
                                                                 pressure = pressure.loc[current_hour],
                                                                 method = 'nrel_numpy',
                                                                 temperature = temp_air.loc[current_hour],
                                                                 )

        # determine the sun apparent elevation corresponding to the azimuth where the horizon profile is provided
        hor_elev = np.interp(sun_positions_rise_set['azimuth'], hor_data.index, hor_data['horizon_elevation'])
        sun_above_horizon = np.array(sun_positions_rise_set['apparent_elevation'] > hor_elev)

        # get the effective sunrise and sunset time, accounting for the horizon
        sun_rise_set_year.loc[day, 'real_sunrise'] = pd.Series(precise_time).where(sun_above_horizon, pd.NaT).dropna().min()  # sunrise time is the earliest time when the sun is above the horizon
        sun_rise_set_year.loc[day, 'real_sunset'] = pd.Series(precise_time).where(sun_above_horizon, pd.NaT).dropna().max()   # sunset time is the last time the sun is above the horizon

    # 7.4 - GET DNI HORIZON LOSS
    dni_loss_hor = pd.Series(index = local_time, dtype=float)

    for hours in local_time:
        hour_start = hours
        hour_end = hour_start + pd.Timedelta(hours=1)
        current_day = hours.normalize()

        if hour_end <= sun_rise_set_year.loc[current_day, 'real_sunrise'] or hour_start >= sun_rise_set_year.loc[current_day, 'real_sunset']:      # if the timestep ends before sunrise or starts after sunset, the DNI is lost
            dni_loss_hor.loc[hours] = 1
        elif hour_start >= sun_rise_set_year.loc[current_day, 'real_sunrise'] and hour_end <= sun_rise_set_year.loc[current_day, 'real_sunset']:   # if the timestamp start and finish are within the sunrise and sunset time, no DNI loss at all
            dni_loss_hor.loc[hours] = 0
        elif hour_start < sun_rise_set_year.loc[current_day, 'real_sunrise'] < hour_end:
            dni_loss_hor.loc[hours] = (sun_rise_set_year.loc[current_day, 'real_sunrise'] - hour_start) / pd.Timedelta(hours=1)
        elif hour_start < sun_rise_set_year.loc[current_day, 'real_sunset'] < hour_end:
            dni_loss_hor.loc[hours] = (hour_end - sun_rise_set_year.loc[current_day, 'real_sunset']) / pd.Timedelta(hours=1)

    dni_hor = dni * (1-dni_loss_hor)

    # 7.5 - DHI HORIZON LOSS
    # Compute the delta across the horizon azimuthal data (it is the 'Base' of the integration shape used later).
    # Note this process will create an array with one value less than hor_elev_angles.
    # Hence a final value is added to restore the removed data and complete the full rotation 0-360 deg
    hor_data['delta_azimuth'] = np.append(np.diff(hor_data.index), 360 - hor_data.index[-1])

    # compute the projected area of the horizon over the sky dome by integrating all the horizon pieces
    # the area of a rectangular sphere sector in steradians from equator is delta_azimuth * (1 - cos(elevation_angle))
    # since the area is not rectangular but trapezoidal, we calculate the area using the right and left elevation, then we average
    area_horizon_strip_rx = np.sum(hor_data['delta_azimuth']/360*2*np.pi * (1-tools.cosd(hor_data['horizon_elevation'])))
    area_horizon_strip_lx = np.sum(hor_data['delta_azimuth']/360*2*np.pi * (1-tools.cosd(np.roll(hor_data['horizon_elevation'], 1))))

    area_horizon_strip = (area_horizon_strip_rx+area_horizon_strip_lx)/2

    # compute the area of the sky dome in steradians (it is the surface of half a sphere)
    area_sky_dome = 2 * np.pi

    # the DHI horizon loss is the ratio between the area of the horizon and the area of the sky dome
    dhi_hor = dhi * (1-area_horizon_strip/area_sky_dome)

    # 7.6 - RECOMPUTE GHI FROM THE COMPONENTS AFTER HORIZON LOSS
    ghi_hor = dhi_hor + dni_hor * tools.cosd(solar_position['zenith'])


    # =================================================================================================
    # 8 - PV FIELD LOOP
    # =================================================================================================
    # Everything that depends on the mounting structure, the module, the inverter and the electrical layout
    # of the field is computed here, once per PV field, from the plane of array irradiance down to the exit
    # of the field MV cabin.

    # hourly AC power delivered by each field at the exit of its MV cabin, in W. Used later to walk the MV tree
    field_power = {}

    # energy at the boundary of each block of the results, per PV field, in Wh. They are the basis of the
    # Performance Ratio, which is computed as a ratio of energies and not as a product of the loss rows
    energy_block_1 = {}   # end of BLOCK 1: exit of the field MV cabin
    energy_block_2 = {}   # end of BLOCK 2: after PV availability and after the MV line leaving the field

    # [kWh/m2] global irradiation on the collector plane of each field, front side, single shed.
    # It is the reference of the Performance Ratio, and the equivalent of the 'GlobInc' of a PVsyst report
    glob_inc = {}

    # [-] optical to electrical conversion coefficient of the rear side of each field, i.e. the share of the
    # rear irradiance the modules actually convert. It is not a loss and does not belong to the loss chain,
    # but it is needed by the consistency check at the end of the run
    bifaciality_factor = {}

    for field in pv_project.pv_fields:

        # =============================================================================================
        # 8.1 - UNPACK FIELD INPUT AND PREPROCESS
        # =============================================================================================

        field_name = field.name

        # POWER
        p_dc_peak = field.power.p_dc_peak * 1000 * 1000      # [W] field peak DC power
        p_ac_nom = field.power.p_ac_nom * 1000 * 1000        # [VA] field nominal inverter power

        # MODULE MOUNTING STRUCTURE
        mount_type = field.mms.mount_type
        mod_configuration = field.mms.mod_configuration
        modules_per_string = field.mms.modules_per_string
        pitch = field.mms.pitch
        tilt_angle = field.mms.tilt_angle
        axis_azimuth = field.mms.axis_azimuth
        backtrack = field.mms.backtrack
        ground_gap = field.mms.ground_gap

        # GROUND
        # if the TMY carries an albedo timeseries, it overrides the albedo declared by the user for this field
        albedo = tmy_albedo if tmy_albedo is not None else field.ground.albedo
        slope_tilt = np.degrees(np.arctan(field.ground.slope_tilt/100))   # convert the ground slope from percent to degrees
        slope_azimuth = field.ground.slope_azimuth

        # MODULE
        Pnom = field.module.Pnom
        mod_height = field.module.height
        mod_width = field.module.width
        cells_in_series = field.module.cells_in_series
        cells_in_parall = field.module.cells_in_parallel
        I_L_ref = field.module.I_L_ref
        imp_ref = field.module.imp_ref
        vmp_ref = field.module.vmp_ref
        alpha_sc = field.module.alpha_sc/1000                # convert from mA/Celsius to A/Celsius
        gamma_ref = field.module.gamma_ref
        mu_gamma = field.module.mu_gamma
        R_sh_ref = field.module.R_sh_ref
        R_sh_0 = field.module.R_sh_0
        R_sh_exp = field.module.R_sh_exp
        R_s = field.module.R_s
        alpha_absorption = field.module.alpha_absorption
        bifaciality = field.module.bifaciality
        substring_type = field.module.substring_type

        # INVERTER
        OND_import = field.inverter.OND_import
        p_ac_inv_nom = field.inverter.Pnom_ac * 1000         # [W] inverter nominal AC output power
        p_ac_inv_max = field.inverter.P_max_ac * 1000        # [W] inverter max AC output power
        v_min_mppt = field.inverter.v_min_mppt
        v_max_mppt = field.inverter.v_max_mppt
        i_max_mppt = field.inverter.i_max_mppt
        n_mppt = field.inverter.n_mppt
        v_out_inv = field.inverter.v_out
        T_inv_p_nom = field.inverter.T_p_nom
        T_inv_p_max = field.inverter.T_p_max

        if OND_import:
            v_eff1 = field.inverter.v_eff1
            v_eff2 = field.inverter.v_eff2
            v_eff3 = field.inverter.v_eff3
            p_in_v_eff1 = np.asarray(field.inverter.p_in_v_eff1) * 1000     # convert from kW to W
            p_in_v_eff2 = np.asarray(field.inverter.p_in_v_eff2) * 1000
            p_in_v_eff3 = np.asarray(field.inverter.p_in_v_eff3) * 1000
            p_out_v_eff1 = np.asarray(field.inverter.p_out_v_eff1) * 1000
            p_out_v_eff2 = np.asarray(field.inverter.p_out_v_eff2) * 1000
            p_out_v_eff3 = np.asarray(field.inverter.p_out_v_eff3) * 1000
            T_eff_1 = field.inverter.T_lim_1
            T_eff_2 = field.inverter.T_lim_2
            p_out_t1 = field.inverter.p_out_T_lim_1 * 1000
            p_out_t2 = field.inverter.p_out_T_lim_2 * 1000
        else:
            inv_euro_eff = field.inverter.inv_euro_eff

        # LV/MV TRAFO
        v_mv = field.mv_trafo.v_out * 1000                   # [V] MV network voltage
        p_trafo_max = field.mv_trafo.p_max * 1000 * 1000     # [VA] max LV/MV transformer power
        trafo_iron_loss = field.mv_trafo.iron_loss
        trafo_copper_loss = field.mv_trafo.copper_loss

        # LOSSES
        iam_model = field.loss.iam_model
        rear_shade_factor = field.loss.rear_shade_factor/100
        thermal_loss_factor = field.loss.thermal_loss_factor
        dc_ohmic_loss_STC = field.loss.dc_ohmic_loss_STC
        ac_lv_cable_len = field.loss.ac_lv_cable_len
        ac_lv_cable_sect = field.loss.ac_lv_cable_sect
        ac_lv_cable_material = field.loss.ac_lv_cable_material
        ac_mv_cable_len = field.loss.ac_mv_cable_len
        ac_mv_cable_sect = field.loss.ac_mv_cable_sect
        ac_mv_cable_material = field.loss.ac_mv_cable_material

        # derived quantities
        n_mod = int(mod_configuration[0])                    # number of modules across the MMS PV table
        stc_eff = Pnom/1000/mod_height/mod_width             # STC efficiency of the module
        is_tracker = 'tracker' in mount_type
        is_portrait = (mod_configuration[1] == 'p')
        is_halfcut = 'half' in substring_type or 'twin' in substring_type

        # total number of modules and total module area of the field. They are used to convert the specific
        # irradiance results [Wh/m2] and the single module results [W] into a field equivalent energy [Wh],
        # so that the losses of different fields can be aggregated into the 'PV project' column
        n_mod_field = p_dc_peak/Pnom
        area_field = n_mod_field * mod_height * mod_width

        if iam_model in ['pan', 'user', 'user_defined']:
            theta_ref = field.module.theta_ref
            iam_ref = field.module.iam_ref


        # =============================================================================================
        # 8.2 - ALLOCATE PV STRUCTURE DATA AND GET ARRAY ORIENTATION
        # =============================================================================================

        # compute the basic geometric data
        if is_portrait:
            pvrow_width = mod_height * n_mod   # width of the PV table when in Portrait
        else:
            pvrow_width = mod_width * n_mod    # width of the PV table when in Landscape

        # for a tracker the tilt angle is a (min, max) tuple, so the largest absolute rotation is taken
        if is_tracker:
            max_tilt = max(abs(angle) for angle in tilt_angle)
        else:
            max_tilt = tilt_angle

        pvrow_height = ground_gap + pvrow_width/2 * tools.sind(max_tilt)   # height of the PV table, measured at the table midline (equivalent to the torque tube axis for a tracker)
        gcr = pvrow_width / pitch

        # Calculate the tracker axis tilt in the global reference frame when on a sloped plane.
        # Axis tilt is the inclination of the rotation axis with respect to horizontal, from 0 deg (horizontal axis) to 90 deg (vertical axis)
        axis_tilt = tracking.calc_axis_tilt(slope_azimuth = slope_azimuth,
                                            slope_tilt = slope_tilt,
                                            axis_azimuth = axis_azimuth,
                                            )

        # Calculate the angle, relative to horizontal, of the line formed by the intersection between the slope
        # containing the tracker axes and a plane perpendicular to the tracker axes.
        # See details here: https://pvlib-python.readthedocs.io/en/stable/reference/generated/pvlib.shading.shaded_fraction1d.html
        cross_axis_tilt = tracking.calc_cross_axis_tilt(slope_azimuth = slope_azimuth,
                                                        slope_tilt = slope_tilt,
                                                        axis_azimuth = axis_azimuth,
                                                        axis_tilt = axis_tilt,
                                                        )

        if is_tracker:
            mount = pvsystem.SingleAxisTrackerMount(axis_tilt = axis_tilt,
                                                    axis_azimuth = axis_azimuth,
                                                    max_angle = tilt_angle,
                                                    backtrack = backtrack,
                                                    gcr = gcr,
                                                    racking_model = 'open_rack',
                                                    module_height = pvrow_height,
                                                    cross_axis_tilt = cross_axis_tilt,
                                                    )
        else:
            # relevant for Fixed tilt only. The azimuth convention is degrees east of north (North = 0, East = 90, South = 180, West = 270)
            surface_azimuth = axis_azimuth + 90
            mount = pvsystem.FixedMount(surface_tilt = tilt_angle,
                                        surface_azimuth = surface_azimuth,
                                        racking_model = 'open_rack',
                                        module_height = pvrow_height,
                                        )

        # compute the structure orientation data (surface tilt and azimuth)
        orientation = mount.get_orientation(solar_zenith = solar_position['apparent_zenith'],
                                            solar_azimuth = solar_position['azimuth'],
                                            )

        # compute the angle of incidence for front and back
        aoi_front = irradiance.aoi(surface_tilt = orientation['surface_tilt'],
                                   surface_azimuth = orientation['surface_azimuth'],
                                   solar_zenith = solar_position['apparent_zenith'],
                                   solar_azimuth = solar_position['azimuth'],
                                   )

        aoi_back = 180 - aoi_front


        # =============================================================================================
        # 8.3 - GET POA IRRADIATION FOR A SINGLE SHED USING THE PEREZ MODEL (SAME AS PVSYST)
        # =============================================================================================

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

        # compute the POA gain over the front side.
        # The resulting irradiation is also stored on its own: it is the reference the Performance Ratio is
        # built on, and it matches the 'GlobInc' of a PVsyst report, i.e. the global irradiation on the
        # collector plane, front side, before the horizon and the near shading losses
        record_loss(Results.POA_GAIN, field_name,
                    np.nansum(ghi) * area_field,
                    np.nansum(irradPOA_single_shed['poa_global']) * area_field)

        glob_inc[field_name] = np.nansum(irradPOA_single_shed['poa_global'])/1000   # [kWh/m2]


        # =============================================================================================
        # 8.4 - POA HORIZON LOSS
        # =============================================================================================
        # the horizon losses on the horizontal components have already been computed once for the project.
        # Here the reduced components are simply transposed onto the plane of array of this field

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

        record_loss('Horizon loss', field_name,
                    np.nansum(irradPOA_single_shed['poa_global']) * area_field,
                    np.nansum(irradPOA_hor['poa_global']) * area_field)


        # =============================================================================================
        # 8.5 - BACK IRRADIANCE
        # =============================================================================================

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


        # =============================================================================================
        # 8.6 - SHADING FACTORS
        # =============================================================================================

        # 8.6.1 - SHADING FACTOR - DIRECT RADIATION - FRONT AND BACK
        # calculate the shaded fraction of the module row width due to self shading
        SF_direct_front = shading.shaded_fraction1d(solar_zenith = solar_position['apparent_zenith'],
                                                    solar_azimuth = solar_position['azimuth'],
                                                    axis_azimuth = axis_azimuth,
                                                    shaded_row_rotation = orientation['surface_tilt'],
                                                    collector_width = pvrow_width,
                                                    pitch = pitch,
                                                    axis_tilt = axis_tilt,
                                                    surface_to_axis_offset = Constants.surface_to_axis_offset,
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
                                                   surface_to_axis_offset = Constants.surface_to_axis_offset,
                                                   cross_axis_slope = cross_axis_tilt,
                                                   shading_row_rotation = 180 - orientation['surface_tilt'],
                                                   )

        # 8.6.2 - SHADING FACTORS - DIFFUSE SKY RADIATION - FRONT AND BACK
        # since we have POA components, to compute the shading loss we need to know the difference of view factor
        # (module -> sky) between the single shed and the multi shed case. The ratio between the two gives the loss
        # of POA irradiance caused by the multi shed layout

        VF_PV_sky_front_multi_row = bifacial.utils.vf_row_sky_2d_integ(surface_tilt = orientation['surface_tilt'],
                                                                       gcr = gcr,
                                                                       x0 = 0,
                                                                       x1 = 1,
                                                                       )

        VF_PV_sky_back_multi_row = bifacial.utils.vf_row_sky_2d_integ(surface_tilt = 180 - orientation['surface_tilt'],
                                                                      gcr = gcr,
                                                                      x0 = 0,
                                                                      x1 = 1,
                                                                      )


        VF_PV_sky_front_single_row = (1 + tools.cosd(orientation['surface_tilt']))/2   # never zero for a tilt below 180 deg

        # View factor for sky is clipped to zero because formula can give slighly negative values for some angles instead of zero.
        # Numerical zero is used instead of zero because view factor is denominator in future fomrmulas. A value equal to zero would resul in cose error
        if is_tracker:
            # surface_tilt is a time series: the floor is applied timestep by timestep
            VF_PV_sky_back_single_row = np.maximum((1 + tools.cosd(180 - orientation['surface_tilt']))/2, Constants.numerical_zero)
        else:
            # surface_tilt is a constant: a single comparison is enough
            VF_PV_sky_back_single_row = max((1 + tools.cosd(180 - orientation['surface_tilt']))/2, Constants.numerical_zero)

        # remove the negative values PVlib returns at the extremes of the tilt range, which are rounding artefacts
        VF_PV_sky_front_multi_row = np.maximum(VF_PV_sky_front_multi_row, 0)
        VF_PV_sky_back_multi_row = np.maximum(VF_PV_sky_back_multi_row, 0)

        # 8.6.3 - SHADING FACTOR - DIFFUSE GROUND RADIATION - FRONT AND BACK
        """
        Since we have POA component, to compute shading loss, we need to know the difference of view factor (module->ground) from the single shed to the multi shed case
        Respect to single shed case, there is also a reduced radiation from sky to ground, since the module will prevent the ground to receive diffuse radiation.
        We will acount this reduction with the View factor Sky-ground for multi shed case VF_gnd_sky_multi_row

        Also, differently from single shed case, the direct-like radiation could also be blocked by the module!
        Hence we need to compute the portion of the ground that receives diffuse radiation only (because shaded by the modules) and the portion receiving direct+circumsolar+diffuse

        Therefore, we also need to know the circumsolar + direct component sum on the horizontal plane
        """

        # compute the number of PV rows to consider in the integration, so that the integration step is less than 5 deg
        max_rows = np.ceil(pvrow_height / (pitch * tools.tand(Constants.max_shading_integration_step)))

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

        # calculate the view factor from module surface to ground, front and back, for the multiple row case
        VF_PV_ground_front_multi_row = bifacial.utils.vf_row_ground_2d_integ(surface_tilt = orientation['surface_tilt'],
                                                                             gcr = gcr,
                                                                             x0 = 0,
                                                                             x1 = 1,
                                                                             )

        VF_PV_ground_back_multi_row = bifacial.utils.vf_row_ground_2d_integ(surface_tilt = 180 - orientation['surface_tilt'],
                                                                            gcr = gcr,
                                                                            x0 = 0,
                                                                            x1 = 1,
                                                                            )

         # View factor for ground is clipped to zero because formula can give slighly negative values for some angles instead of zero.
         # Numerical zero is used instead of zero because view factor is denominator in future fomrmulas. A value equal to zero would resul in cose error
        if is_tracker:
            # surface_tilt is a time series: the floor is applied timestep by timestep
            VF_PV_ground_front_single_row = np.maximum((1 - tools.cosd(orientation['surface_tilt']))/2, Constants.numerical_zero)
        else:
            # surface_tilt is a constant: a single comparison is enough
            VF_PV_ground_front_single_row = max((1 - tools.cosd(orientation['surface_tilt']))/2, Constants.numerical_zero)

        # the back side looks at (1+cos(tilt))/2, which is never zero for a tilt below 180 deg, hence no floor needed
        VF_PV_ground_back_single_row = (1 - tools.cosd(180 - orientation['surface_tilt']))/2

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

        # compute the circumsolar and direct component on the horizontal plane, assuming a flat surface
        irrad_horizontal = get_irradiance_full_components(surface_tilt = 0,
                                                          surface_azimuth = orientation['surface_azimuth'],
                                                          solar_zenith = solar_position['apparent_zenith'],
                                                          solar_azimuth = solar_position['azimuth'],
                                                          aoi = solar_position['apparent_zenith'],   # for a flat surface, the incidence angle equals the sun zenith angle
                                                          ghi = ghi_hor,
                                                          dhi = dhi_hor,
                                                          dni = dni_hor,
                                                          dni_extra = dni_extra,
                                                          airmass = airmass_relative,
                                                          albedo = 0,   # not needed here, since only the ratio between direct+circumsolar and diffuse is used
                                                          )

        direct_fraction = (irrad_horizontal['poa_direct'] + irrad_horizontal['poa_sky_diffuse_circumsolar'])/irrad_horizontal['poa_global']
        diffuse_fraction = 1 - direct_fraction


        # =============================================================================================
        # 8.6.4 - GROUND BRIGHTNESS SEEN BY THE MODULES
        # =============================================================================================
        
        # The ground has very non uniform brightness; in fact some areas (e.g. under the module) are shadowed, while other areas (between rows) are lit by sun.
        # By calculating the average of ground brightness, I would overestimate the ground component.
        # Hence, ground component is calcualted by brightness weighted averaged of each ground section with its brightness. View factor and ground brighness are caluclated for each ground strip and then integrated all together.
        # However, for tracker the view factor an brightness chances at each angle. therefore, we calculate such data for 15 standard angle. All other tracker position are calcualted via interpolation of such 7 standard positions
        # This is the same method used by PVsyst

        # Compute the standard tilt angles on which the geometry is pre-computed. A tracker sweeps its whole
        # range, for a fixed structure a single angle is enough
        if is_tracker:
            standard_angles = np.linspace(0, np.nanmax(orientation['surface_tilt']), Constants.ground_tilt_steps)
        else:
            standard_angles = np.array([float(np.nanmax(orientation['surface_tilt']))])

        # Pre-compute, for every standard angle at once, the weights of the view of the ground from each face.
        # The ground positions are the same for every angle, and front and rear come out of a single call
        x_ground, weights_ground_front, weights_ground_back = vf_row_ground_2d_weights(standard_angles = standard_angles,
                                                                                       gcr = gcr,
                                                                                       pvrow_height = pvrow_height,
                                                                                       pitch = pitch,
                                                                                       )

        # Pre-compute the sky seen from every ground point. It does not depend on the face, only on the structure.
        # PVlib returns the result shaped (ground points, angles), so it is transposed to match the weights above
        vf_sky_ground = np.asarray(utils.vf_ground_sky_2d(rotation = standard_angles,
                                                          gcr = gcr,
                                                          x = x_ground,
                                                          pitch = pitch,
                                                          height = pvrow_height,
                                                          max_rows = Constants.ground_max_rows_sky,
                                                          )).T

        # interpolate the pre-computed quantities on the actual tilt of every timestep.
        # For a fixed structure the tilt is a constant, for a tracker it is a time series: np.broadcast_to
        # turns both into an array of the right length without duplicating any memory
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
        tan_psi = np.asarray(utils._solar_projection_tangent(solar_zenith = solar_position['apparent_zenith'],
                                                             solar_azimuth = solar_position['azimuth'],
                                                             surface_azimuth = orientation['surface_azimuth'],
                                                             ), dtype=float)
        tan_psi = np.nan_to_num(tan_psi, nan=0.0, posinf=0.0, neginf=0.0)

        # THE SHADOW OF THE ROWS ON THE GROUND, resolved at every timestep.
        # Each edge of the row is projected onto the ground along the direction of the sun: a point at height h
        # lands at h*tan_psi further along the ground. With the right edge of the row higher, the front looks to
        # the left, so a sun on the front side casts the shadow towards the right.
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


        # =============================================================================================
        # 8.7 - NEAR SHADING LOSSES
        # =============================================================================================

        # create the structure to store irradPOA_near_shade (front + back)
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

        # POA direct radiation after near shading, front and back
        irradPOA_near_shade['poa_front_direct'] = irradPOA_hor['poa_direct'] * (1-SF_direct_front)
        irradPOA_near_shade['poa_back_direct'] = irradPOA_back['poa_direct'] * (1-SF_direct_back) * (1-rear_shade_factor)

        # POA sky isotropic diffuse radiation after near shading, front and back
        irradPOA_near_shade['poa_front_sky_diffuse_isotropic'] = irradPOA_hor['poa_sky_diffuse_isotropic'] / VF_PV_sky_front_single_row * VF_PV_sky_front_multi_row
        irradPOA_near_shade['poa_back_sky_diffuse_isotropic'] = irradPOA_back['poa_sky_diffuse_isotropic'] / VF_PV_sky_back_single_row * VF_PV_sky_back_multi_row * (1-rear_shade_factor)

        # POA sky circumsolar diffuse radiation after near shading. Circumsolar is assumed to behave like direct radiation, so the shading factor is applied
        irradPOA_near_shade['poa_front_sky_diffuse_circumsolar'] = irradPOA_hor['poa_sky_diffuse_circumsolar'] * (1-SF_direct_front)
        irradPOA_near_shade['poa_back_sky_diffuse_circumsolar'] = irradPOA_back['poa_sky_diffuse_circumsolar'] * (1-SF_direct_back) * (1-rear_shade_factor)

        # POA sky horizon diffuse radiation after near shading
        # this is set to zero assuming that modules are not able to see horizon strip. In fact horizon strip bad is 6.5° according to Perez. Such band is likely hidden by next PV structure
        # in previous caluclation, it was assumed that horizon view factor was the same as sky view factor. However, such assumption is completely wrong and produced overestimated horizon irradiance
        irradPOA_near_shade['poa_front_sky_diffuse_horizon'] = 0  # irradPOA_near_shade['poa_front_sky_diffuse_horizon'] = irradPOA_hor['poa_sky_diffuse_horizon'] / VF_PV_sky_front_single_row * VF_PV_sky_front_multi_row
        irradPOA_near_shade['poa_back_sky_diffuse_horizon'] = 0  # irradPOA_near_shade['poa_back_sky_diffuse_horizon'] = irradPOA_back['poa_sky_diffuse_horizon'] / VF_PV_sky_back_single_row * VF_PV_sky_back_multi_row * (1-rear_shade_factor)
        

        # compute the summary of the diffuse components
        irradPOA_near_shade['poa_front_sky_diffuse'] = irradPOA_near_shade['poa_front_sky_diffuse_isotropic'] + irradPOA_near_shade['poa_front_sky_diffuse_circumsolar'] + irradPOA_near_shade['poa_front_sky_diffuse_horizon']
        irradPOA_near_shade['poa_back_sky_diffuse'] = irradPOA_near_shade['poa_back_sky_diffuse_isotropic'] + irradPOA_near_shade['poa_back_sky_diffuse_circumsolar'] + irradPOA_near_shade['poa_back_sky_diffuse_horizon']

        # POA ground diffuse radiation after near shading, front and back.
        # The chain is: (irradiance leaving the open ground) x (ground brightness reduction) x (view factor row -> visible ground)
        #   - dividing by the single row view factor gives back the irradiance leaving an unobstructed ground (ghi*albedo)
        #   - the multi row view factor accounts for how much of that ground the module actually sees
        #   - ground_brightness front and back is an averaged reduction factor that trasforms irradiance leaving unobstructed ground (ghi*albedo) to the actual irradiance leaving the ground under the module
        
        irradPOA_near_shade['poa_front_ground_diffuse'] = irradPOA_hor['poa_ground_diffuse'] / VF_PV_ground_front_single_row * VF_PV_ground_front_multi_row * ground_brightness_front
        irradPOA_near_shade['poa_back_ground_diffuse'] = irradPOA_back['poa_ground_diffuse'] / VF_PV_ground_back_single_row * VF_PV_ground_back_multi_row * ground_brightness_back * (1-rear_shade_factor)

        # PREVIOUS FORMULATION - DO NOT DELETE
        # Previously, instead of using the averaged ground_brightness, it was used (diffuse_fraction * VF_gnd_sky_multi_row + direct_fraction * fraction_gnd_direct)
        #   - diffuse_fraction: fraction of diffuse radiation in POA radiation (=1-direct_fraction)
        #   - VF_gnd_sky_multi_row: View factor from ground to sky
        #   - direct_fraction: fraction of direct radiation in POA radiation (=1-diffuse_fraction)
        #   - fraction_gnd_direct: fraction of area under PV structures that receives direct radiation
        # irradPOA_near_shade['poa_front_ground_diffuse'] = irradPOA_hor['poa_ground_diffuse'] / VF_PV_ground_front_single_row * VF_PV_ground_front_multi_row * (diffuse_fraction * VF_gnd_sky_multi_row + direct_fraction * fraction_gnd_direct)
        # irradPOA_near_shade['poa_back_ground_diffuse'] = irradPOA_back['poa_ground_diffuse'] / VF_PV_ground_back_single_row * VF_PV_ground_back_multi_row * (diffuse_fraction * VF_gnd_sky_multi_row + direct_fraction * fraction_gnd_direct) * (1-rear_shade_factor)

        # complete the irradPOA_near_shade variable with the remaining components
        irradPOA_near_shade['poa_front_diffuse'] = irradPOA_near_shade['poa_front_sky_diffuse'] + irradPOA_near_shade['poa_front_ground_diffuse']
        irradPOA_near_shade['poa_back_diffuse'] = irradPOA_near_shade['poa_back_sky_diffuse'] + irradPOA_near_shade['poa_back_ground_diffuse']
        irradPOA_near_shade['poa_front'] = irradPOA_near_shade['poa_front_direct'] + irradPOA_near_shade['poa_front_diffuse']
        irradPOA_near_shade['poa_back'] = irradPOA_near_shade['poa_back_direct'] + irradPOA_near_shade['poa_back_diffuse']
        irradPOA_near_shade['poa_global'] = irradPOA_near_shade['poa_front'] + irradPOA_near_shade['poa_back']

        # To respect the PVsyst convention, the shading loss is computed respect to the front side only.
        # However, irradPOA_near_shade already contains the components for the back side including the shading loss. Those components will be used to calculate the POA gain on the back side later
        record_loss('Near shading loss', field_name,
                    np.nansum(irradPOA_hor['poa_global']) * area_field,
                    np.nansum(irradPOA_near_shade['poa_front']) * area_field)


        # =============================================================================================
        # 8.8 - IAM LOSS
        # =============================================================================================

        # 8.8.1 - DIRECT COMPONENT IAM
        # if the IAM model is GlassAR, the IAM angles are generated and the relative IAM values are computed
        if iam_model in ['glassar', 'glass', 'glass ar', 'ar glass', 'glass_ar', 'ar_glass']:
            theta_ref = np.arange(0, 90, 0.5)
            iam_ref = iam.martin_ruiz(aoi = theta_ref,
                                      a_r = Constants.iam_a_r,
                                      )

        # the IAM direct component loss is calculated by interpolating the angle of incidence over the IAM loss function
        iam_direct_front = iam.interp(aoi = aoi_front,
                                      theta_ref = theta_ref,
                                      iam_ref = iam_ref,
                                      method = 'cubic',
                                      normalize = True,
                                      )

        iam_direct_back = iam.interp(aoi = aoi_back,
                                     theta_ref = theta_ref,
                                     iam_ref = iam_ref,
                                     method = 'cubic',
                                     normalize = True,
                                     )

        # 8.8.2 - DIFFUSE COMPONENT IAM
        # compute the IAM loss function for the diffuse radiation
        iam_function = interp1d(theta_ref, iam_ref, kind='cubic', bounds_error=False, fill_value=(iam_ref[0], iam_ref[-1]))

        # compute the IAM loss for each separate diffuse component.
        iam_diffuse_front_sky = iam.marion_integrate(function = iam_function,
                                                     surface_tilt = np.max(orientation['surface_tilt']),
                                                     region = 'sky',
                                                     num = None,
                                                     )

        iam_diffuse_front_horizon = iam.marion_integrate(function = iam_function,
                                                         surface_tilt = np.max(orientation['surface_tilt']),
                                                         region = 'horizon',
                                                         num = None,
                                                         )

        iam_diffuse_front_ground = iam.marion_integrate(function = iam_function,
                                                        surface_tilt = np.max(orientation['surface_tilt']),
                                                        region = 'ground',
                                                        num = None,
                                                        )

        iam_diffuse_back_sky = iam.marion_integrate(function = iam_function,
                                                    surface_tilt = 180 - np.max(orientation['surface_tilt']),
                                                    region = 'sky',
                                                    num = None,
                                                    )

        iam_diffuse_back_horizon = iam.marion_integrate(function = iam_function,
                                                        surface_tilt = 180 - np.max(orientation['surface_tilt']),
                                                        region = 'horizon',
                                                        num = None,
                                                        )

        iam_diffuse_back_ground = iam.marion_integrate(function = iam_function,
                                                       surface_tilt = 180 - np.max(orientation['surface_tilt']),
                                                       region = 'ground',
                                                       num = None,
                                                       )

        # preallocate a new DataFrame for irradPOA_iam, then calculate the components after the IAM loss
        irradPOA_iam = pd.DataFrame(columns = POA_columns)

        # IAM loss for the direct radiation
        irradPOA_iam['poa_front_direct'] = irradPOA_near_shade['poa_front_direct'] * iam_direct_front
        irradPOA_iam['poa_back_direct'] = irradPOA_near_shade['poa_back_direct'] * iam_direct_back

        # IAM for the sky diffuse radiation components. For circumsolar, the IAM loss of the direct component is applied
        irradPOA_iam['poa_front_sky_diffuse_isotropic'] = irradPOA_near_shade['poa_front_sky_diffuse_isotropic'] * iam_diffuse_front_sky
        irradPOA_iam['poa_back_sky_diffuse_isotropic'] = irradPOA_near_shade['poa_back_sky_diffuse_isotropic'] * iam_diffuse_back_sky
        irradPOA_iam['poa_front_sky_diffuse_circumsolar'] = irradPOA_near_shade['poa_front_sky_diffuse_circumsolar'] * iam_direct_front
        irradPOA_iam['poa_back_sky_diffuse_circumsolar'] = irradPOA_near_shade['poa_back_sky_diffuse_circumsolar'] * iam_direct_back
        irradPOA_iam['poa_front_sky_diffuse_horizon'] = irradPOA_near_shade['poa_front_sky_diffuse_horizon'] * iam_diffuse_front_horizon
        irradPOA_iam['poa_back_sky_diffuse_horizon'] = irradPOA_near_shade['poa_back_sky_diffuse_horizon'] * iam_diffuse_back_horizon

        # compute the summary of the sky diffuse components
        irradPOA_iam['poa_front_sky_diffuse'] = irradPOA_iam['poa_front_sky_diffuse_isotropic'] + irradPOA_iam['poa_front_sky_diffuse_circumsolar'] + irradPOA_iam['poa_front_sky_diffuse_horizon']
        irradPOA_iam['poa_back_sky_diffuse'] = irradPOA_iam['poa_back_sky_diffuse_isotropic'] + irradPOA_iam['poa_back_sky_diffuse_circumsolar'] + irradPOA_iam['poa_back_sky_diffuse_horizon']

        # IAM for the ground diffuse radiation components
        irradPOA_iam['poa_front_ground_diffuse'] = irradPOA_near_shade['poa_front_ground_diffuse'] * iam_diffuse_front_ground
        irradPOA_iam['poa_back_ground_diffuse'] = irradPOA_near_shade['poa_back_ground_diffuse'] * iam_diffuse_back_ground

        # complete the irradPOA_iam variable with the remaining components
        irradPOA_iam['poa_front_diffuse'] = irradPOA_iam['poa_front_sky_diffuse'] + irradPOA_iam['poa_front_ground_diffuse']
        irradPOA_iam['poa_back_diffuse'] = irradPOA_iam['poa_back_sky_diffuse'] + irradPOA_iam['poa_back_ground_diffuse']
        irradPOA_iam['poa_front'] = irradPOA_iam['poa_front_direct'] + irradPOA_iam['poa_front_diffuse']
        irradPOA_iam['poa_back'] = irradPOA_iam['poa_back_direct'] + irradPOA_iam['poa_back_diffuse']
        irradPOA_iam['poa_global'] = irradPOA_iam['poa_front'] + irradPOA_iam['poa_back']

        # To respect the PVsyst convention, the IAM loss is computed respect to the front side only.
        # However, irradPOA_iam already contains the components for the back side including the shading loss and IAM. Those components will be used to calculate the POA gain on the back side later
        record_loss('IAM loss', field_name,
                    np.nansum(irradPOA_near_shade['poa_front']) * area_field,
                    np.nansum(irradPOA_iam['poa_front']) * area_field)


        # =============================================================================================
        # 8.9 - SOILING LOSS
        # =============================================================================================
        # the effective irradiance is the irradiance hitting the PV module after IAM, soiling and snow loss.
        # It is the radiation that can actually be converted into energy. Snow loss should be added here in the future

        soiling_factor = 1 + RESULTS.loc['Soiling loss', field_name]

        # multiply every item of irradPOA_iam by the soiling loss
        irradPOA_soiling = soiling_factor * irradPOA_iam

        # To respect the PVsyst convention, the Soiling loss is computed respect to the front side only.
        # However, irradPOA_soiling already contains the components for the back side including the shading loss, IAM and soiling. Those components will be used to calculate the POA gain on the back side later
        record_loss('Soiling loss', field_name,
                    np.nansum(irradPOA_iam['poa_front']) * area_field,
                    np.nansum(irradPOA_soiling['poa_front']) * area_field)
        
        
        # =============================================================================================
        # 8.10 - SNOW LOSSES - TO BE IMPLEMENTED IN THE FUTURE
        # =============================================================================================
        """
        if latitude <=55:
            snow_loss = 0     # soiling loss expresses in % - default is 2
        else:
            snow_loss = max(3, (latitude-55)/5*3)
        """

        # =============================================================================================
        # 8.11 - POA GAIN BACK SIDE, BIFACIALITY EFFICIENCY AND EFFECTIVE IRRADIANCE
        # =============================================================================================
        # apply the bifaciality to the back irradiance. Every front component of irradPOA_effective is the same as
        # irradPOA_soiling. Bifaciality is the ratio between the STC efficiency of the front side and of the back side

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
        irradPOA_effective['poa_back'] = irradPOA_soiling['poa_back'] * bifaciality

        irradPOA_effective['poa_global'] = irradPOA_effective['poa_front'] + irradPOA_effective['poa_back']

        # share of the total plane of array irradiance the modules can convert, once the rear side has been
        # de-rated by the bifaciality factor. Together with the STC efficiency it is what turns the optical
        # chain into the electrical one, exactly as the 'PV conversion' step of a PVsyst loss diagram
        bifaciality_factor[field_name] = (np.nansum(irradPOA_effective['poa_global'])
                                          / np.nansum(irradPOA_soiling['poa_front'] + irradPOA_soiling['poa_back']))


        # NOTE: this is the OPTICAL gain of the rear side, before the bifaciality factor, so that it can be
        # compared one to one with the 'Global Irradiance on rear side' line of a PVsyst loss diagram.
        # The bifaciality factor is an optical-to-electrical conversion coefficient, applied above together with
        # the STC efficiency, exactly as PVsyst does: it is therefore NOT part of the loss chain
        record_loss('POA gain back', field_name,
                    np.nansum(irradPOA_soiling['poa_front']) * area_field,
                    np.nansum(irradPOA_soiling['poa_front'] + irradPOA_soiling['poa_back']) * area_field)


        # =============================================================================================
        # 8.12 - IRRADIANCE AND TEMPERATURE LOSS - PVSYST METHOD
        # =============================================================================================
        # this section computes the loss due to operating at non-STC conditions, for both temperature and
        # irradiance, following the same approach used by PVsyst

        # 8.12.1 - REAL MODULE CONDITIONS
        # first the cell temperature is computed using the PVsyst method
        temp_cell = temperature.pvsyst_cell(poa_global = irradPOA_effective['poa_global'],
                                            temp_air = temp_air,
                                            wind_speed = wind_speed,
                                            u_c = thermal_loss_factor,                                # combined heat loss factor coefficient. The default value of 29 is representative of freestanding modules with the rear surface exposed to open air (e.g. rack mounted)
                                            u_v = 0.0,                                                # combined heat loss factor influenced by wind
                                            module_efficiency = Constants.module_efficiency_thermal,  # module external efficiency. It would be Pdc/(poa_global * module_area), however Pdc is not known a priori. The PVsyst default value is 0.1
                                            alpha_absorption = alpha_absorption,                      # PV module absorption coefficient = 1 - reflectivity
                                            )

        # then the parameters of the single diode equation are calculated
        pvsyst_params = pvsystem.calcparams_pvsyst(effective_irradiance = irradPOA_effective['poa_global'],
                                                   temp_cell = temp_cell,
                                                   alpha_sc = alpha_sc,
                                                   gamma_ref = gamma_ref,
                                                   mu_gamma = mu_gamma,
                                                   I_L_ref = I_L_ref,
                                                   I_o_ref = Constants.I_o_ref,
                                                   R_sh_ref = R_sh_ref,
                                                   R_sh_0 = R_sh_0,
                                                   R_s = R_s,
                                                   cells_in_series = cells_in_series,
                                                   R_sh_exp = R_sh_exp,
                                                   EgRef = Constants.EgRef,
                                                   irrad_ref = Constants.irrad_ref,
                                                   temp_ref = Constants.temp_ref,
                                                   )

        # the function above returns the results in a tuple. They are rearranged into a DataFrame for convenience
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
                                                   method = 'chandrupatla',   # either 'chandrupatla', 'brentq' or 'newton'. 'chandrupatla' is slow and guaranteed to converge, Brent's method is mid-speed and guaranteed to converge (not working here), 'newton' is fast but NOT guaranteed to converge
                                                   )

        module_real_power = max_power_point['p_mp']

        # 8.12.2 - NON IDEALITY LOSS CALCULATION
        # the theoretical power is calculated as STC efficiency * effective irradiance
        #
        # NOTE: this is the point where the loss chain switches from irradiance [Wh/m2] to electrical energy [Wh],
        # and the conversion factor is stc_eff. It is therefore the ONLY row of the chain whose input does not
        # numerically equal the output of the previous row, and that is by design, not an inconsistency:
        # the ratio between the two is exactly stc_eff. It used to carry the bifaciality factor as well, which
        # has now been moved to its own row at the end of section 8.10
        theoretical_power = irradPOA_effective['poa_global'] * mod_height * mod_width * stc_eff

        record_loss('Temp + Irr loss', field_name,
                    np.nansum(theoretical_power) * n_mod_field,
                    np.nansum(module_real_power) * n_mod_field)


        # =============================================================================================
        # 8.13 - PARTIAL SHADING - ELECTRICAL EFFECTS LOSS - EMPIRIC CALCULATION
        # =============================================================================================
        # Electrical effects happen when a substring is partially shaded, causing irradiation mismatch and hotspots.
        # Substrings that are fully shaded have already been accounted for in the shading loss.
        # Here we account for those substrings that are partially shaded and partially unshaded.
        # In this method, the partial shading loss is calculated proportionally to the shaded area of the module

        # calculate the number of substrings and the width [m] of a module cell (in the direction of pvrow_width)
        if is_portrait and is_halfcut:          # Portrait oriented half-cut module
            total_substrings = cells_in_parall*n_mod
            cell_width = mod_height/(cells_in_series/Constants.substrings_per_module)

        elif is_portrait and not is_halfcut:    # Portrait oriented full-cut module
            total_substrings = n_mod
            cell_width = mod_height/(cells_in_series/Constants.substrings_per_module)

        else:                                   # Landscape oriented module, both half-cut and full-cut
            total_substrings = Constants.substrings_per_module*n_mod
            cell_width = mod_width/(2*Constants.substrings_per_module)

        # calculate the shaded fraction of the partially shaded substring (only the decimal part of the total shaded area)
        extra_shaded_substring = SF_direct_front*total_substrings - np.floor(SF_direct_front*total_substrings)

        # calculate the size of one cell relative to the size of one substring
        cell_relative_width = cell_width/(pvrow_width/total_substrings)

        # calculate the power loss due to the electrical effect:
        # if the shaded part is shorter than the length of a cell, the substring power loss is proportional to the shaded area of the cell
        # if the shaded part is longer than the length of a cell, the whole direct component of that substring is lost
        power_loss_partial_shading = np.where(extra_shaded_substring < cell_relative_width,
                                              extra_shaded_substring/cell_relative_width/total_substrings*(irradPOA_effective['poa_front_direct'] + irradPOA_effective['poa_front_sky_diffuse_circumsolar'])/irradPOA_effective['poa_global'],
                                              1/total_substrings*(irradPOA_effective['poa_front_direct'] + irradPOA_effective['poa_front_sky_diffuse_circumsolar'])/irradPOA_effective['poa_global']
                                              )

        record_loss('Partial shading loss', field_name,
                    np.nansum(module_real_power) * n_mod_field,
                    np.nansum(module_real_power * (1-power_loss_partial_shading)) * n_mod_field)

        # compute the total loss inside the PV module (electrical effect, module quality, LID, mismatch front, mismatch back)
        module_internal_loss = ((1-power_loss_partial_shading)
                                * (1 + RESULTS.loc['Module quality gain', field_name])
                                * (1 + RESULTS.loc['LID + LeTID loss', field_name])
                                * (1 + RESULTS.loc['Mismatch front loss', field_name])
                                * (1 + RESULTS.loc['Mismatch back loss', field_name])
                                )

        # the static module losses are chained onto the running energy, so that they can also be aggregated
        # into the 'PV project' column as an energy weighted average
        e_run = np.nansum(module_real_power * (1-power_loss_partial_shading)) * n_mod_field
        for row in ['Module quality gain', 'LID + LeTID loss', 'Mismatch front loss', 'Mismatch back loss']:
            e_next = e_run * (1 + RESULTS.loc[row, field_name])
            record_loss(row, field_name, e_run, e_next)
            e_run = e_next


        # =============================================================================================
        # 8.14 - STRINGING
        # =============================================================================================
        """
        PVyield uses a basic stringing mechanism that consist of 2xinverter configuration, one 'big' inverter with more PV strings and one 'small' inverter with less PV string.
        An inverter configuration is defined by:
            - the number of inverter of that configuration
            - the number of string connected to that inverter
            - the MPPT configuration for that inverter
        All inverter configuration infos above are included in the dictionary InverterConfig

        The number of inverter with more PV strings and less PV strings are calculated to match the total number of string in the field.
        E.g. if the field has 107 PV strings to split among 12 inverter, the simulation will assume:
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

        the code will result in:
        n_mppt_less_string = 4   n_less_string = 2
        n_mppt_more_string = 8   n_more_string = 3

        NOTE: n_more_input is always n_less_input + 1
        """

        n_inv_field = int(np.ceil(p_ac_nom/p_ac_inv_nom))       # number of inverters of the field. p_ac_nom was already rounded to an integer number of inverters during the input check
        tot_string_field = p_dc_peak/modules_per_string/Pnom    # total number of PV strings in the field

        # n_string_base is the number of strings every inverter carries, n_extra_string is how many inverters carry one more
        n_string_base, n_extra_string = divmod(tot_string_field, n_inv_field)

        InverterConfig = {'big': {},
                          'small': {},
                          }

        if n_extra_string == 0:
            # the strings divide evenly: every inverter carries the same number of strings and there is
            # no 'small' configuration at all. Its number of inverters is set to zero, so that it carries no weight
            InverterConfig['big']['n_string'] = n_string_base
            InverterConfig['big']['n_inv'] = n_inv_field
            InverterConfig['small']['n_string'] = n_string_base
            InverterConfig['small']['n_inv'] = 0
        else:
            # n_extra_string inverters carry one string more than the others
            InverterConfig['big']['n_string'] = n_string_base + 1
            InverterConfig['big']['n_inv'] = n_extra_string
            InverterConfig['small']['n_string'] = n_string_base
            InverterConfig['small']['n_inv'] = n_inv_field - n_extra_string

        # safety net: the strings distributed among the inverters must add up to the strings of the field.
        # Should a future change break the distribution again, this stops the run instead of quietly losing energy
        distributed_strings = (InverterConfig['big']['n_string'] * InverterConfig['big']['n_inv']
                               + InverterConfig['small']['n_string'] * InverterConfig['small']['n_inv'])
        if distributed_strings != tot_string_field:
            raise RuntimeError(f'PV field "{field_name}": internal error while distributing the strings among the inverters. '
                               f'{distributed_strings} strings were allocated instead of {tot_string_field}')

        # distribute the strings of each inverter among its MPPT inputs, following the same logic
        for config in ['big', 'small']:
            n_string_config = InverterConfig[config]['n_string']

            # n_less_string is the number of strings every MPPT input carries, n_mppt_more_string is how many carry one more
            n_less_string, n_mppt_more_string = divmod(int(n_string_config), int(n_mppt))

            InverterConfig[config]['n_less_string'] = n_less_string                      # number of strings connected to the MPPT with less strings connected
            InverterConfig[config]['n_more_string'] = n_less_string + 1                  # number of strings connected to the MPPT with more strings connected
            InverterConfig[config]['n_mppt_more_string'] = n_mppt_more_string            # number of MPPT inputs with the max amount of strings connected
            InverterConfig[config]['n_mppt_less_string'] = n_mppt - n_mppt_more_string   # number of MPPT inputs with the minimum amount of strings connected

            # same safety net at MPPT level
            distributed_mppt = (InverterConfig[config]['n_more_string'] * InverterConfig[config]['n_mppt_more_string']
                                + InverterConfig[config]['n_less_string'] * InverterConfig[config]['n_mppt_less_string'])
            if distributed_mppt != n_string_config:
                raise RuntimeError(f'PV field "{field_name}": internal error while distributing the strings among the MPPT inputs '
                                   f'of the "{config}" inverter. {distributed_mppt} strings were allocated instead of {n_string_config}')


        # =============================================================================================
        # 8.15 - INVERTER MPPT FINDING
        # =============================================================================================

        # 8.15.1 - PREALLOCATE VARIABLES
        inverter_columns = ['vmp_string',             # [V] string voltage at module max power point conditions (module MPP voltage * modules per string), before applying the MPPT limitations
                            'imp_string',             # [A] module current at module max power point conditions, before applying the MPPT limitations
                            'pmp_module',             # [W] module power at module max power point conditions, before applying the MPPT limitations
                            'i_max_mppt',             # [A] max MPPT current limit applicable to the string (MPPT current limit divided by the number of strings in that MPPT)
                            'vmp_mppt_i_max_limit',   # [V] MPPT input lower voltage threshold due to the max current limit of each MPPT input
                            'vmp_mppt_i_inv_limit',   # [V] MPPT input lower voltage threshold as per inverter specification
                            'vmp_mppt_low_limit',     # [V] actual MPPT input lower voltage threshold, the maximum between vmp_mppt_i_max_limit and vmp_mppt_i_inv_limit
                            'vmp_mppt_high_limit',    # [V] MPPT input upper voltage threshold as per inverter specification
                            'vmp_mppt',               # [V] MPPT input voltage, after all the voltage constraints have been accounted for
                            'imp_mppt_ideal',         # [A] MPPT DC current of each string, before accounting for the internal PV losses (LID, quality loss, mismatch front, mismatch back, partial shading)
                            'imp_mppt',               # [A] MPPT DC current of each string, after accounting for the internal PV losses
                            'pmp_mppt_module',        # [W] power of the PV module at inverter MPPT conditions
                            'dc_ohm_powr_loss',       # [Wh] DC ohmic loss in the strings
                            'dc_ohm_loss',            # [-] ohmic loss in the string, in the [0-1] range
                            'p_dc_string_out',        # [W] DC power at the terminals of the PV modules, before the string cable loss
                            'p_dc_inv_in',            # [W] inverter DC input power, after the string cable loss
                            'p_ac_inv_out',           # [W] inverter AC output power
                            'p_dc_inv_in_ideal',      # [W] inverter DC input power for an unconstrained inverter
                            ]

        for config in ['big', 'small']:
            InverterConfig[config]['mppt_less_input'] = pd.DataFrame(index = local_time, columns = inverter_columns)
            InverterConfig[config]['mppt_more_input'] = pd.DataFrame(index = local_time, columns = inverter_columns)

            for mppt in ['mppt_less_input', 'mppt_more_input']:
                InverterConfig[config][mppt]['vmp_string'] = max_power_point['v_mp'] * modules_per_string
                InverterConfig[config][mppt]['imp_string'] = max_power_point['i_mp'] * module_internal_loss
                InverterConfig[config][mppt]['pmp_module'] = max_power_point['p_mp'] * module_internal_loss
                InverterConfig[config][mppt]['vmp_mppt_i_inv_limit'] = v_min_mppt
                InverterConfig[config][mppt]['vmp_mppt_high_limit'] = v_max_mppt

            # DC power the strings would deliver at the TRUE maximum power point of the modules, i.e. with an
            # MPPT free of any voltage or current constraint. It is taken at the module terminals, before the
            # string cable, and it is the reference of the 'Inverter MPPT V/I loss' row further below
            InverterConfig[config]['p_dc_inv_in_ideal'] = InverterConfig[config]['mppt_more_input']['pmp_module'] * modules_per_string * InverterConfig[config]['n_string']

        # 8.15.2 - COMPUTE MPPT VOLTAGE BOUNDARIES
        """
        The inverter MPPT input has a maximum current limit. This limit is for the whole MPPT.
        So if there are multiple strings connected to that MPPT input, each string can carry a maximum value equal to the inverter MPPT current limit divided by the strings connected to the MPPT.
        e.g. if the inverter MPPT input has a limit of 30A, and can connect 3 strings max, each string can carry 10A max

        We need to be sure that the MPP of the module generates a current that is lower than this limit (10A in the example above).
        Therefore, the inverter MPPT voltage scan cannot select a voltage corresponding to a current higher than this limit.
        Since the IV curve is monotonically decreasing, this max current limit generates a lower voltage limit for the inverter MPPT scan.

        In other words, the inverter MPPT max current limit generates a MPPT voltage lower limit. If the voltage of the string is lower than this limit, the current produced will be higher than the acceptable inverter MPPT current limit
        The voltage limit due to current is compared with the voltage limit of the MPPT, and the most stringent of the two (so the maximum among the two) is taken into account as MPPT voltage scan limit
        """

        for config in ['big', 'small']:
            for mppt, n_string_key in [('mppt_less_input', 'n_less_string'), ('mppt_more_input', 'n_more_string')]:

                if InverterConfig[config][n_string_key] == 0:
                    # no string connected to this MPPT input type, hence no current limit to compute
                    InverterConfig[config][mppt]['vmp_mppt_i_max_limit'] = 0
                    continue

                # starting from the current limit of the inverter MPPT, compute the max current limit applicable to a single string
                InverterConfig[config][mppt]['i_max_mppt'] = i_max_mppt/module_internal_loss/InverterConfig[config][n_string_key]

                # if the max current limit exceeds the photocurrent of the module, the limit is set to the photocurrent itself.
                # In that case the v_from_i function returns a voltage equal to 0, i.e. no constraint
                InverterConfig[config][mppt]['i_max_mppt'] = np.where(InverterConfig[config][mppt]['i_max_mppt'] >= iv_params['photocurrent'],
                                                                      iv_params['photocurrent'],
                                                                      InverterConfig[config][mppt]['i_max_mppt'])

                # compute the lower boundary of the MPPT voltage due to the max current limit
                InverterConfig[config][mppt]['vmp_mppt_i_max_limit'] = modules_per_string * pvsystem.v_from_i(current = InverterConfig[config][mppt]['i_max_mppt'],
                                                                                                             photocurrent = iv_params['photocurrent'],
                                                                                                             saturation_current = iv_params['saturation_current'],
                                                                                                             resistance_series = iv_params['resistance_series'],
                                                                                                             resistance_shunt = iv_params['resistance_shunt'],
                                                                                                             nNsVth = iv_params['nNsVth'],
                                                                                                             method = 'lambertw',   # 'lambertw', 'chandrupatla' or 'brentq'. NOTE: 'brentq' is limited to the 1st quadrant only
                                                                                                             )

        # 8.15.3 - RECOMPUTE THE MPPT VOLTAGE TO MAKE SURE IT IS INSIDE THE BOUNDARIES
        for config in ['big', 'small']:
            for mppt in ['mppt_less_input', 'mppt_more_input']:

                # the actual lower voltage limit is the most stringent between the current driven limit and the inverter specification
                InverterConfig[config][mppt]['vmp_mppt_low_limit'] = np.where(InverterConfig[config][mppt]['vmp_mppt_i_max_limit'] > v_min_mppt,
                                                                              InverterConfig[config][mppt]['vmp_mppt_i_max_limit'],
                                                                              v_min_mppt)

                # compute the actual MPP voltage by clipping the original string voltage to the lower and upper MPPT voltage boundaries
                InverterConfig[config][mppt]['vmp_mppt'] = np.clip(InverterConfig[config][mppt]['vmp_string'],
                                                                   InverterConfig[config][mppt]['vmp_mppt_low_limit'],
                                                                   InverterConfig[config][mppt]['vmp_mppt_high_limit'])

                # 8.14.4 - COMPUTE THE MPPT CURRENT FROM THE VOLTAGE
                InverterConfig[config][mppt]['imp_mppt_ideal'] = pvsystem.i_from_v(voltage = InverterConfig[config][mppt]['vmp_mppt']/modules_per_string,
                                                                                   photocurrent = iv_params['photocurrent'],
                                                                                   saturation_current = iv_params['saturation_current'],
                                                                                   resistance_series = iv_params['resistance_series'],
                                                                                   resistance_shunt = iv_params['resistance_shunt'],
                                                                                   nNsVth = iv_params['nNsVth'],
                                                                                   method = 'lambertw',   # 'lambertw', 'newton' or 'brentq'. NOTE: 'brentq' is limited to the 1st quadrant only
                                                                                   )

                # the actual current flowing in the string is lower than the result of the IV curve, because of the module
                # internal losses (LID, quality loss, front mismatch, back mismatch, partial shading)
                InverterConfig[config][mppt]['imp_mppt'] = InverterConfig[config][mppt]['imp_mppt_ideal'] * module_internal_loss

                # the PV module power at MPPT conditions is calculated
                InverterConfig[config][mppt]['pmp_mppt_module'] = InverterConfig[config][mppt]['imp_mppt'] * InverterConfig[config][mppt]['vmp_mppt'] / modules_per_string


        # =============================================================================================
        # 8.16 - DC OHMIC LOSS
        # =============================================================================================

        for config in ['big', 'small']:
            for suffix, n_string_key, n_mppt_key, mppt in [('less', 'n_less_string', 'n_mppt_less_string', 'mppt_less_input'),
                                                           ('more', 'n_more_string', 'n_mppt_more_string', 'mppt_more_input')]:

                # compute the equivalent resistance needed to compute the DC string cable losses
                InverterConfig[config]['eq_resistance_'+suffix+'_input'] = pvsystem.dc_ohms_from_percent(vmp_ref = vmp_ref,                     # voltage at maximum power in reference conditions [V]
                                                                                                        imp_ref = imp_ref,                     # current at maximum power in reference conditions [A]
                                                                                                        dc_ohmic_percent = dc_ohmic_loss_STC,
                                                                                                        modules_per_string = modules_per_string,
                                                                                                        strings = InverterConfig[config][n_string_key],   # number of parallel strings in the array
                                                                                                        )

                if InverterConfig[config][n_string_key] == 0:
                    # no strings connected to this MPPT input type, hence no DC losses
                    InverterConfig[config][mppt]['dc_ohm_powr_loss'] = 0
                else:
                    # current flowing through one MPPT input: the string current times the strings connected to it
                    i_mppt_input = InverterConfig[config][mppt]['imp_mppt'] * InverterConfig[config][n_string_key]

                    InverterConfig[config][mppt]['dc_ohm_powr_loss'] = pvsystem.dc_ohmic_losses(resistance = InverterConfig[config]['eq_resistance_'+suffix+'_input'],
                                                                                                current = i_mppt_input,
                                                                                                )

                    # scale the loss of a single MPPT input up to all the MPPT inputs of that type
                    InverterConfig[config][mppt]['dc_ohm_powr_loss'] = InverterConfig[config][mppt]['dc_ohm_powr_loss'] * InverterConfig[config][n_mppt_key]

                # power leaving the PV modules connected to this MPPT input type, at the module terminals
                InverterConfig[config][mppt]['p_dc_string_out'] = InverterConfig[config][mppt]['pmp_mppt_module'] * modules_per_string * InverterConfig[config][n_string_key] * InverterConfig[config][n_mppt_key]


                # NOTE ON THE OPERATING POINT: the MPPT search above does not 'feel' the cable resistance, therefore maximises the power at the MODULE terminals and ignores the cable resistance.
                # A real inverter would 'feel' the cable resistance, finding MPPT points at its own terminaly. Therefore, a real inverter MPPT would be at slightly higher voltage and drawing a slightly lower current.
                # This effect is deliberately left out: it is a second order effect accounting for a fraction (<0.01 %) of yield
                # It is possible to add this effect by adding the cable resistance to R_s before calcparams_pvsyst, which also changes the cell temperature and the MPPT search, and is not worth it for a gain of this size.
                InverterConfig[config][mppt]['p_dc_inv_in'] = (InverterConfig[config][mppt]['p_dc_string_out'] - InverterConfig[config][mppt]['dc_ohm_powr_loss'])

            InverterConfig[config]['dc_ohm_powr_loss'] = InverterConfig[config]['mppt_less_input']['dc_ohm_powr_loss'] + InverterConfig[config]['mppt_more_input']['dc_ohm_powr_loss']
            InverterConfig[config]['p_dc_string_out'] = InverterConfig[config]['mppt_less_input']['p_dc_string_out'] + InverterConfig[config]['mppt_more_input']['p_dc_string_out']
            InverterConfig[config]['p_dc_inv_in'] = InverterConfig[config]['mppt_less_input']['p_dc_inv_in'] + InverterConfig[config]['mppt_more_input']['p_dc_inv_in']

        # Total power at the module terminals and total DC power entering the inverters of the field.
        p_dc_string_out = InverterConfig['big']['p_dc_string_out'] * InverterConfig['big']['n_inv'] + InverterConfig['small']['p_dc_string_out'] * InverterConfig['small']['n_inv']

        # MPPT VOLTAGE AND CURRENT CONSTRAINTS
        # The MPPT of the inverter cannot always sit on the maximum power point of the modules: its voltage
        # window and its maximum input current push the operating point away from it. The energy lost that way
        # is the difference between the power the strings deliver at the constrained operating point and the
        # power they would deliver at the true maximum power point, both taken at the module terminals.
        #
        # This loss used to be hidden inside a row called 'Inverter thresholds loss', computed as the
        # ratio between the IDEAL clipped output and the REAL pre-clipping one. Mixing the two chains gave a
        # row that did not measure its own stage, and with a real OND file it even came out POSITIVE: on the
        # reference project it reported +0.500%, i.e. a gain, where the inverter can only take energy away.
        # The reason is that the ideal chain starts from the power at the module terminals while the real one
        # starts after the string cable, so the row was silently reporting the DC cable loss with the sign
        # reversed. The two effects are now separated and each row measures exactly the stage it names:
        #   'Inverter MPPT V/I loss'  : here, on the DC side, where the constraint physically acts
        #   'Inverter clipping loss'  : further below, on the AC side, where the output is capped
        # Both are always negative or zero, and the loss chain telescopes without needing any correction.
        p_dc_string_out_ideal = InverterConfig['big']['p_dc_inv_in_ideal'] * InverterConfig['big']['n_inv'] + InverterConfig['small']['p_dc_inv_in_ideal'] * InverterConfig['small']['n_inv']

        record_loss('Inverter MPPT V/I loss', field_name,
                    np.nansum(p_dc_string_out_ideal),
                    np.nansum(p_dc_string_out))
        p_dc_inv_in = InverterConfig['big']['p_dc_inv_in'] * InverterConfig['big']['n_inv'] + InverterConfig['small']['p_dc_inv_in'] * InverterConfig['small']['n_inv']

        record_loss('DC ohmic wiring loss', field_name,
                    np.nansum(p_dc_string_out),
                    np.nansum(p_dc_inv_in))


        # =============================================================================================
        # 8.17 - INVERTER EFFICIENCY
        # =============================================================================================
        # compute the output power as a function of the three efficiency curves specified in the OND:
        # P_out = f(P_in, V_in, efficiency curves)

        if OND_import:
            for config in ['big', 'small']:

                # average MPPT voltage of the inverter, weighted on the number of MPPT inputs of each type
                v_dc_in = (InverterConfig[config]['mppt_more_input']['vmp_mppt'] * InverterConfig[config]['n_mppt_more_string']
                           + InverterConfig[config]['mppt_less_input']['vmp_mppt'] * InverterConfig[config]['n_mppt_less_string']) / n_mppt

                InverterConfig[config]['p_ac_inv_out'] = inv_Pout_from_Pin(v_eff1 = v_eff1,
                                                                           v_eff2 = v_eff2,
                                                                           v_eff3 = v_eff3,
                                                                           p_in_v_eff1 = p_in_v_eff1,
                                                                           p_in_v_eff2 = p_in_v_eff2,
                                                                           p_in_v_eff3 = p_in_v_eff3,
                                                                           p_out_v_eff1 = p_out_v_eff1,
                                                                           p_out_v_eff2 = p_out_v_eff2,
                                                                           p_out_v_eff3 = p_out_v_eff3,
                                                                           p_dc_inv_in = InverterConfig[config]['p_dc_inv_in'],
                                                                           v_dc_in = v_dc_in,
                                                                           )

            p_ac_inv_out_preclipping = InverterConfig['big']['p_ac_inv_out'] * InverterConfig['big']['n_inv'] + InverterConfig['small']['p_ac_inv_out'] * InverterConfig['small']['n_inv']

        else:
            # with no OND file, a flat european efficiency is used
            p_ac_inv_out_preclipping = p_dc_inv_in * inv_euro_eff/100

        record_loss('Inverter efficiency loss', field_name,
                    np.nansum(p_dc_inv_in),
                    np.nansum(p_ac_inv_out_preclipping))


        # =============================================================================================
        # 8.18 - INVERTER CLIPPING LOSS
        # =============================================================================================
        # calculate the max inverter output power as a function of temperature. The inverter temperature is
        # assumed to be equal to the air temperature (outdoor installation)

        if OND_import:
            inv_temp_ref = [T_inv_p_max, T_inv_p_nom, T_eff_1, T_eff_2]
            inv_p_ref_t_corrected = [p_ac_inv_max, p_ac_inv_nom, p_out_t1, p_out_t2]
        else:
            inv_temp_ref = [T_inv_p_max, T_inv_p_nom]
            inv_p_ref_t_corrected = [p_ac_inv_max, p_ac_inv_nom]

        # the reference points are sorted by increasing temperature, so that the interpolation works regardless
        # of the order in which the temperature thresholds are declared in the OND or in the Inverter class
        sorted_points = sorted(zip(inv_temp_ref, inv_p_ref_t_corrected))
        inv_temp_ref = [point[0] for point in sorted_points]
        inv_p_ref_t_corrected = [point[1] for point in sorted_points]
        fill_value = (inv_p_ref_t_corrected[0], inv_p_ref_t_corrected[-1])

        # create the power = f(temperature) interpolation function
        inv_max_p_out_t_function = interp1d(inv_temp_ref, inv_p_ref_t_corrected, kind='slinear', bounds_error=False, fill_value=fill_value)

        # apply the interpolation function to the air temperature to find the max inverter power of the whole field
        inv_max_p_out_t_corrected = pd.Series(inv_max_p_out_t_function(temp_air), index=local_time) * n_inv_field

        p_ac_inv_out = p_ac_inv_out_preclipping.clip(upper = inv_max_p_out_t_corrected)

        # NOTE: as in the previous version, the threshold loss compares the IDEAL clipped output against the REAL
        # pre-clipping output, so that it also embeds the MPPT constraint losses
        # the clipping is measured on the REAL chain alone, so the row is always negative or zero
        record_loss('Inverter clipping loss', field_name,
                    np.nansum(p_ac_inv_out_preclipping),
                    np.nansum(p_ac_inv_out))


        # =============================================================================================
        # 8.19 - CURTAILMENT AT INVERTER OUTPUT
        # =============================================================================================
        # first of the two possible curtailment points at PV field level. The AC power of the whole field is
        # limited right at the inverter output, before the LV AC cable

        p_field = p_ac_inv_out

        if field.field_curtailment.grid_curtail_point == 'inverter out':
            curtail_power = field.field_curtailment.grid_curtail_power * 1000 * 1000   # [VA]
            p_field_curtailed = p_field.clip(upper = curtail_power)
            record_loss('Curtailment @ inverter', field_name, np.nansum(p_field), np.nansum(p_field_curtailed))
            p_field = p_field_curtailed
        else:
            record_loss('Curtailment @ inverter', field_name, np.nansum(p_field), np.nansum(p_field))


        # =============================================================================================
        # 8.20 - AC LV OHMIC LOSS
        # =============================================================================================
        # compute the losses in the cabling from each inverter to its LV/MV transformer

        p_in_ac_lv_cable = p_field/n_inv_field   # power flowing through a single LV cable

        p_ac_lv_loss = ac_ohm_loss(p_ac = p_in_ac_lv_cable,
                                   v_p2p = v_out_inv,
                                   leng = ac_lv_cable_len,
                                   material = ac_lv_cable_material,
                                   section = ac_lv_cable_sect,
                                   )

        p_out_ac_lv_cable = p_in_ac_lv_cable * (1-p_ac_lv_loss)

        record_loss('LV AC cable loss', field_name,
                    np.nansum(p_in_ac_lv_cable) * n_inv_field,
                    np.nansum(p_out_ac_lv_cable) * n_inv_field)

        p_field = p_out_ac_lv_cable * n_inv_field


        # =============================================================================================
        # 8.21 - LV/MV TRAFO LOSS
        # =============================================================================================

        n_trafo = np.ceil(p_ac_nom/p_trafo_max)   # number of LV/MV transformers in the field
        p_trafo = p_ac_nom/n_trafo                # actual rating of each transformer

        p_trafo_in = p_field / n_trafo            # power entering a single transformer

        p_trafo_out = transformer.simple_efficiency(input_power = p_trafo_in,
                                                    no_load_loss = trafo_iron_loss/100,
                                                    load_loss = trafo_copper_loss/100,
                                                    transformer_rating = p_trafo,
                                                    )

        record_loss('LV/MV trafo loss', field_name,
                    np.nansum(p_trafo_in) * n_trafo,
                    np.nansum(p_trafo_out) * n_trafo)

        p_field = p_trafo_out * n_trafo


        # =============================================================================================
        # 8.22 - AUXILIARIES LOSS
        # =============================================================================================
        # the auxiliaries (trackers, SCADA, lighting, security, HVAC of the cabins) are supplied from the MV side
        # of the LV/MV transformer, hence their consumption is subtracted between the transformer and the MV cable

        e_in = np.nansum(p_field)
        p_field = p_field * (1 + RESULTS.loc['Auxiliaries loss', field_name])
        record_loss('Auxiliaries loss', field_name, e_in, np.nansum(p_field))


        # =============================================================================================
        # 8.23 - AC MV OHMIC LOSS
        # =============================================================================================
        # compute the losses in the MV cabling connecting the transformer antennas to the PV field collection cabin

        p_in_ac_mv_cable = p_field / n_trafo   # power flowing through a single MV antenna

        # The auxiliaries loss is now applied to p_field beforehand, so the formula is unambiguous
        p_ac_mv_loss = ac_ohm_loss(p_ac = p_in_ac_mv_cable,
                                   v_p2p = v_mv,
                                   leng = ac_mv_cable_len,
                                   material = ac_mv_cable_material,
                                   section = ac_mv_cable_sect,
                                   )

        p_out_ac_mv_cable = p_in_ac_mv_cable * (1-p_ac_mv_loss)

        record_loss('MV AC cable loss', field_name,
                    np.nansum(p_in_ac_mv_cable) * n_trafo,
                    np.nansum(p_out_ac_mv_cable) * n_trafo)

        p_field = p_out_ac_mv_cable * n_trafo


        # =============================================================================================
        # 8.24 - FIRST YEAR DEGRADATION
        # =============================================================================================
        # the degradation is applied to the hourly power BEFORE the curtailment, so that the curtailment loss
        # is computed on the power the plant will actually deliver during its first year of operation

        e_in = np.nansum(p_field)
        p_field = p_field * (1 + RESULTS.loc['First year degradation', field_name])
        record_loss('First year degradation', field_name, e_in, np.nansum(p_field))


        # =============================================================================================
        # 8.25 - CURTAILMENT AT FIELD MV CABIN
        # =============================================================================================
        # second of the two possible curtailment points at PV field level. The curtailment is applied at the MV
        # cabin of the field, after the MV antennas and before the MV connection line leaving the field

        if field.field_curtailment.grid_curtail_point == 'field mv cabin':
            curtail_power = field.field_curtailment.grid_curtail_power * 1000 * 1000   # [VA]
            p_field_curtailed = p_field.clip(upper = curtail_power)
            record_loss('Curtailment @ MV cabin', field_name, np.nansum(p_field), np.nansum(p_field_curtailed))
            p_field = p_field_curtailed
        else:
            record_loss('Curtailment @ MV cabin', field_name, np.nansum(p_field), np.nansum(p_field))


        # =============================================================================================
        # 8.26 - PV AVAILABILITY
        # =============================================================================================
        # the availability is applied here so that the power walking the MV tree already accounts for it

        # energy leaving the PV field at its MV cabin. It closes BLOCK 1 of the results, and is used further
        # below to compute the Performance Ratio the same way PVsyst does, as a ratio of energies
        energy_block_1[field_name] = np.nansum(p_field)

        e_in = np.nansum(p_field)
        p_field = p_field * (1 + RESULTS.loc['PV availability', field_name])
        record_loss('PV availability', field_name, e_in, np.nansum(p_field))

        # store the hourly AC power leaving this field, used to walk the MV connection tree
        field_power[field_name] = p_field


    # =================================================================================================
    # 9 - MV CONNECTION TREE
    # =================================================================================================
    # The PV fields are linked to each other through MV connection lines, forming a tree whose single exit node
    # is the Point of Connection. The power of a field therefore travels through its own MV line and then through
    # the MV lines of every field standing between it and the exit of the tree.
    #
    # The tree is walked from the exit node backwards: the power entering the line of a field is the sum of the
    # power produced by that field plus the power delivered by every field connected upstream of it.
    #
    # ATTRIBUTION CONVENTION: the 'MV connection loss (field-to-field)' reported in the column of a field is the
    # loss of ITS OWN line only, computed on the whole power flowing through it (own production plus upstream
    # fields). The losses of the lines standing downstream, closer to the Point of Connection, are NOT charged
    # to that field. The value shown in the 'PV project' column, on the contrary, accounts for every line of the tree.

    exit_field = pv_project.exit_field                       # the single field feeding the Point of Connection
    v_mv_project = pv_project.pv_fields[0].mv_trafo.v_out * 1000   # [V] MV voltage, common to all the fields

    # project curtailment settings
    project_curtail_point = pv_project.project_curtailment.grid_curtail_point
    if project_curtail_point is not None:
        project_curtail_power = pv_project.project_curtailment.grid_curtail_power * 1000 * 1000   # [VA]

    line_in = {}     # hourly power entering the MV line of each field, before any curtailment
    line_out = {}    # hourly power leaving the MV line of each field

    def mv_line_loss(pv_field, p_in):
        """
        Applies the ohmic loss of the MV connection line leaving one PV field.

        Args:
            pv_field (PV_Field):  the PV field owning the line
            p_in (pd.Series):     [VA] hourly power entering the line

        Returns:
            p_out (pd.Series):    [VA] hourly power leaving the line
        """

        mv_connection = pv_field.mv_connection

        # no line declared: the two nodes are physically adjacent, hence no loss
        if mv_connection.line_len is None or mv_connection.line_len == 0:
            return p_in

        if mv_connection.auto_sizing:
            raise NotImplementedError(f'PV field "{pv_field.name}": automatic sizing of the MV connection line is not implemented yet. '
                                      f'Please set auto_sizing = False and declare line_nr and line_sect manually')

        # the power is split among the parallel cables of the line
        p_per_cable = p_in / mv_connection.line_nr

        line_loss = ac_ohm_loss(p_ac = p_per_cable,
                                v_p2p = v_mv_project,
                                leng = mv_connection.line_len * 1000,   # convert from km to m
                                material = mv_connection.material,
                                section = mv_connection.line_sect,
                                )

        return p_per_cable * (1-line_loss) * mv_connection.line_nr

    def walk_mv_tree(pv_field):
        """
        Recursively computes the hourly power leaving the MV line of one PV field, walking the tree upstream.

        Args:
            pv_field (PV_Field):  the PV field to be processed

        Returns:
            p_out (pd.Series):    [VA] hourly power leaving the MV line of the field
        """

        if pv_field.name in line_out:
            return line_out[pv_field.name]

        # the power entering the line is the field production plus everything coming from upstream
        p_in = field_power[pv_field.name].copy()
        for upstream_field in pv_field.incoming_fields:
            p_in = p_in + walk_mv_tree(upstream_field)

        line_in[pv_field.name] = p_in

        # the project curtailment 'Final MV line in' applies at the entry of the LAST line of the tree only
        p_line_in = p_in
        if pv_field is exit_field and project_curtail_point == 'final mv line in':
            p_line_in = p_in.clip(upper = project_curtail_power)

        p_out = mv_line_loss(pv_field, p_line_in)

        # the loss reported for this field is the loss of its own line, computed after any curtailment
        record_loss('MV connection loss (field-to-field)', pv_field.name, np.nansum(p_line_in), np.nansum(p_out))

        # energy of THIS field alone at the end of BLOCK 2, i.e. after its availability and after its own MV line.
        # The loss of the line is shared with the fields upstream, so it is applied here as a factor
        line_factor = np.nansum(p_out)/np.nansum(p_line_in) if np.nansum(p_line_in) else 0
        energy_block_2[pv_field.name] = np.nansum(field_power[pv_field.name]) * line_factor

        line_out[pv_field.name] = p_out
        return p_out

    # walk the tree starting from the exit node. Every other field is reached recursively
    p_poc = walk_mv_tree(exit_field)

    # ---------------------------------------------------------------------------------------------
    # aggregate the MV tree results into the 'PV project' column
    # ---------------------------------------------------------------------------------------------
    # The curtailment sits in the middle of the line cascade, so the two effects cannot be separated simply by
    # ordering the rows. They are instead grouped into two exact multiplicative factors whose product equals
    # the overall ratio between the energy delivered at the exit of the tree and the energy produced by the fields.

    e_fields_total = sum(np.nansum(field_power[f.name]) for f in pv_project.pv_fields)   # energy produced by all the fields, at their MV cabins
    e_final_line_in = np.nansum(line_in[exit_field.name])                                # energy entering the last line of the tree, before curtailment
    e_poc = np.nansum(p_poc)                                                             # energy leaving the last line of the tree

    if project_curtail_point == 'final mv line in':
        e_final_line_curtailed = np.nansum(line_in[exit_field.name].clip(upper = project_curtail_power))
    else:
        e_final_line_curtailed = e_final_line_in

    # factor 1: the ohmic losses of EVERY line of the tree, the last one included
    mv_tree_factor = (e_final_line_in/e_fields_total) * (e_poc/e_final_line_curtailed)

    # factor 2: the project curtailment applied at the entry of the last line
    curtail_factor = e_final_line_curtailed/e_final_line_in

    RESULTS.loc['MV connection loss (field-to-field)', Results.PROJECT_COLUMN] = mv_tree_factor - 1
    RESULTS.loc['Curtailment @ Final MV line in', Results.PROJECT_COLUMN] = curtail_factor - 1

    # these two rows are handled explicitly, so they are removed from the automatic energy aggregation
    energy_in.pop('MV connection loss (field-to-field)', None)
    energy_out.pop('MV connection loss (field-to-field)', None)


    # =================================================================================================
    # 10 - COMMON HV PARTS
    # =================================================================================================
    # MV/HV transformer, HV connection line, project curtailment and grid availability. They belong to the shared
    # part of the plant, therefore they are computed once on the aggregated power of all the fields and are
    # reported in the 'PV project' column only.

    def apply_project_curtailment(row, point, p_in):
        """
        Applies the project curtailment at one of the available application points, and stores the resulting loss.

        Args:
            row (str):           name of the RESULTS row where the loss is stored
            point (str):         canonical name of the curtailment application point handled by this call
            p_in (pd.Series):    [VA] hourly power entering the curtailment point

        Returns:
            p_out (pd.Series):   [VA] hourly power after the curtailment
        """

        if project_curtail_point != point:
            RESULTS.loc[row, Results.PROJECT_COLUMN] = 0
            return p_in

        p_out = p_in.clip(upper = project_curtail_power)
        RESULTS.loc[row, Results.PROJECT_COLUMN] = np.nansum(p_out)/np.nansum(p_in) - 1
        return p_out

    # 10.1 - CURTAILMENT AT THE HV TRANSFORMER INPUT (equivalent to the exit of the last MV line)
    p_poc = apply_project_curtailment('Curtailment @ HV trafo in', 'hv trafo in', p_poc)

    # 10.2 - MV/HV TRANSFORMER LOSS
    if pv_project.hv_trafo.p_nom is not None:
        p_hv_trafo_in = p_poc
        p_hv_trafo_out = transformer.simple_efficiency(input_power = p_hv_trafo_in,
                                                       no_load_loss = pv_project.hv_trafo.iron_loss/100,
                                                       load_loss = pv_project.hv_trafo.copper_loss/100,
                                                       transformer_rating = pv_project.hv_trafo.p_nom * 1000 * 1000,
                                                       )

        RESULTS.loc['HV trafo loss', Results.PROJECT_COLUMN] = np.nansum(p_hv_trafo_out)/np.nansum(p_hv_trafo_in) - 1
        p_poc = p_hv_trafo_out
    else:
        # no MV/HV transformer declared: the Point of Connection is on the MV network
        RESULTS.loc['HV trafo loss', Results.PROJECT_COLUMN] = 0

    # 10.3 - CURTAILMENT AT THE HV TRANSFORMER OUTPUT (equivalent to the entry of the HV line)
    p_poc = apply_project_curtailment('Curtailment @ HV trafo out', 'hv trafo out', p_poc)

    # 10.4 - HV CONNECTION LINE LOSS
    # the HV line loss is a flat percentage per km of line, as declared in the HV_Connection class
    if pv_project.hv_connection.line_len is not None:
        hv_line_loss = pv_project.hv_connection.line_len * Constants.hv_line_loss_per_km/100
        RESULTS.loc['HV connection loss', Results.PROJECT_COLUMN] = -hv_line_loss
        p_poc = p_poc * (1-hv_line_loss)
    else:
        RESULTS.loc['HV connection loss', Results.PROJECT_COLUMN] = 0

    # 10.5 - CURTAILMENT AT THE HV LINE OUTPUT
    p_poc = apply_project_curtailment('Curtailment @ HV line out', 'hv line out', p_poc)

    # 10.6 - GRID AVAILABILITY
    # the value has already been stored by Results.assign_results, it is applied here to the hourly power
    p_poc = p_poc * (1 + RESULTS.loc['Grid availability', Results.PROJECT_COLUMN])

    # hourly net energy delivered by the project at the Point of Connection.
    # The simulation timestep is one hour, hence the average power over the step is also the energy of the step
    energy_hourly = p_poc/1000
    energy_hourly.name = 'Net energy at PoC [kWh]'


    # =================================================================================================
    # 11 - BUILD THE PROJECT COLUMN OF THE LOSS CHAIN
    # =================================================================================================
    # Every loss belonging to the PV field block is aggregated across the fields as the ratio between the sum of
    # the energies leaving that stage and the sum of the energies entering it. This is an energy weighted average
    # of the individual field losses, and it keeps the product of the chain consistent with the project Net Yield.

    for row in energy_in:
        if energy_in[row]:
            RESULTS.loc[row, Results.PROJECT_COLUMN] = energy_out[row]/energy_in[row] - 1
        else:
            RESULTS.loc[row, Results.PROJECT_COLUMN] = 0


    # =================================================================================================
    # 12 - COMPUTE PERFORMANCE RATIO AND NET YIELD
    # =================================================================================================
    # The Performance Ratio is a ratio of energies, exactly as PVsyst computes it, and NOT the product of the loss rows
    # See the note in Results.compute_performance for the reason.

    # energy delivered by the whole project at the end of each block, in Wh
    energy_project = {1: sum(energy_block_1.values()),
                      2: sum(energy_block_2.values()),
                      3: np.nansum(energy_hourly)*1000,
                      }

    p_dc_project = sum(f.power.p_dc_peak for f in pv_project.pv_fields) * 1000 * 1000   # [Wp]

    for field in pv_project.pv_fields:
        RESULTS = Results.compute_performance(RESULTS = RESULTS,
                                              column = field.name,
                                              glob_inc = glob_inc[field.name],
                                              energy = {1: energy_block_1[field.name],
                                                        2: energy_block_2[field.name],
                                                        },
                                              p_dc_peak = field.power.p_dc_peak * 1000 * 1000,
                                              )

    # the reference irradiation of the project is the average of the fields, weighted on their peak DC power,
    # so that the project Performance Ratio stays consistent with the Net Yield of the whole plant
    glob_inc_project = (sum(glob_inc[f.name]*f.power.p_dc_peak for f in pv_project.pv_fields)
                        / sum(f.power.p_dc_peak for f in pv_project.pv_fields))

    RESULTS = Results.compute_performance(RESULTS = RESULTS,
                                          column = Results.PROJECT_COLUMN,
                                          glob_inc = glob_inc_project,
                                          energy = energy_project,
                                          p_dc_peak = p_dc_project,
                                          )

    # CONSISTENCY CHECK
    # The Performance Ratio is now a ratio of energies and the loss rows are only descriptive, so the two are
    # no longer forced to agree. They are still tied by an exact relation, which is worth verifying:
    #
    #     PR = product of the loss rows x bifaciality factor
    #
    # because the bifaciality is the only coefficient of the optical to electrical conversion that is not
    # already accounted for by the normalisation on the peak DC power. A mismatch means that a stage of the
    # chain does not start where the previous one ended, which is exactly the kind of error that stays
    # invisible in the results table while quietly moving the numbers
    for field in pv_project.pv_fields:
        chain_product = 1.0
        for row in Results.PV_ROWS:
            value = RESULTS.loc[row, field.name]
            chain_product = chain_product * (1 + (0 if pd.isna(value) else value))

        expected = chain_product * bifaciality_factor[field.name]
        computed = RESULTS.loc[Results.PR_PV_ONLY, field.name]/100

        if computed and abs(expected/computed - 1) > 0.005:
            warnings.warn(f'PV field "{field.name}": the Performance Ratio computed from the delivered energy '
                          f'({computed*100:.2f} %) differs by more than 0.5% from the product of the loss rows '
                          f'corrected by the bifaciality factor ({expected*100:.2f} %). This points to a stage of '
                          f'the loss chain that does not start where the previous one ends, and the results '
                          f'should be reviewed')

    # the same relation is verified on the project column, over the WHOLE chain. The check above only covers
    # the rows of the first block, so on its own it would stay blind to an inconsistency introduced by the MV
    # tree, by the HV parts or by one of the project curtailment points
    chain_product = 1.0
    for row in Results.PV_ROWS + Results.FIELDS_ROWS + Results.PROJECT_ROWS:
        value = RESULTS.loc[row, Results.PROJECT_COLUMN]
        chain_product = chain_product * (1 + (0 if pd.isna(value) else value))

    # the bifaciality of the project is the energy weighted average of the fields, consistently with the way
    # every other row of the project column is aggregated
    energy_weight = {f.name: energy_block_1[f.name] for f in pv_project.pv_fields}
    total_weight = sum(energy_weight.values())
    bifaciality_project = (sum(bifaciality_factor[f.name]*energy_weight[f.name] for f in pv_project.pv_fields)/total_weight
                           if total_weight else 1.0)

    expected = chain_product * bifaciality_project
    computed = RESULTS.loc[Results.PR_PROJECT, Results.PROJECT_COLUMN]/100

    if computed and abs(expected/computed - 1) > 0.005:
        warnings.warn(f'PV project: the Performance Ratio computed from the delivered energy '
                      f'({computed*100:.2f} %) differs by more than 0.5% from the product of the whole loss '
                      f'chain corrected by the bifaciality factor ({expected*100:.2f} %). This points to a '
                      f'stage of the loss chain that does not start where the previous one ends, and the '
                      f'results should be reviewed')


    # =================================================================================================
    # 13 - FORMAT AND PRINT THE RESULTS
    # =================================================================================================

    RESULTS = Results.format_results(RESULTS)

    # print the loss chain of the project on the console. The full table, holding one column per PV field,
    # is written to the Excel report by excel_generator()
    Results.print_results(RESULTS)

    return RESULTS, energy_hourly



# ==========================================================================
# SUPPORT FUNCTIONS
# ==========================================================================
    """#!!! """ # placeholder for support function

def get_irradiance_full_components(surface_tilt,
                                   surface_azimuth,
                                   solar_zenith,
                                   solar_azimuth,
                                   aoi,
                                   ghi,
                                   dhi,
                                   dni,
                                   dni_extra,
                                   airmass,
                                   albedo,
                                   ):
    """
    Computes the full set of irradiance components on a plane of array, using the Perez-Driesse transposition model
    (the same family of models used by PVsyst), returning every component separately instead of the global value only.

    Splitting the diffuse irradiance into its isotropic, circumsolar and horizon parts is required because each part
    behaves differently under near shading (view factors) and under incidence angle losses.

    Args:
        surface_tilt (float or pd.Series):     [deg] tilt of the plane of array, 0 = horizontal facing up, 90 = vertical
        surface_azimuth (float or pd.Series):  [deg] azimuth of the plane of array, measured clockwise from north
        solar_zenith (pd.Series):              [deg] apparent solar zenith angle
        solar_azimuth (pd.Series):             [deg] solar azimuth angle
        aoi (pd.Series):                       [deg] angle of incidence between the sun and the plane of array
        ghi (pd.Series):                       [W/m2] global horizontal irradiance
        dhi (pd.Series):                       [W/m2] diffuse horizontal irradiance
        dni (pd.Series):                       [W/m2] direct normal irradiance
        dni_extra (pd.Series):                 [W/m2] extraterrestrial normal irradiance
        airmass (pd.Series):                   [-] relative airmass
        albedo (float or pd.Series):           [-] ground surface albedo

    Returns:
        components (pd.DataFrame):  DataFrame holding the following columns, all in [W/m2]:
            'poa_global'                  global irradiance on the plane of array
            'poa_direct'                  direct (beam) irradiance on the plane of array
            'poa_diffuse'                 total diffuse irradiance (sky diffuse + ground diffuse)
            'poa_ground_diffuse'          diffuse irradiance reflected by the ground
            'poa_sky_diffuse'             total sky diffuse irradiance (isotropic + circumsolar + horizon)
            'poa_sky_diffuse_isotropic'   isotropic diffuse irradiance, coming from the whole sky dome
            'poa_sky_diffuse_circumsolar' circumsolar diffuse irradiance, coming from the area around the sun
            'poa_sky_diffuse_horizon'     horizon diffuse irradiance, coming from the strip just above the horizon
    """

    # compute the full sky diffuse components (isotropic, circumsolar, horizon)
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

    # the horizon diffuse contribution can become negative in the model. Negative values have no physical meaning and are zeroed
    diffuse_sky['poa_horizon'] = np.where(diffuse_sky['poa_horizon'] < 0, 0, diffuse_sky['poa_horizon'])

    # compute the ground reflected diffuse component
    diffuse_ground = irradiance.get_ground_diffuse(surface_tilt = surface_tilt,
                                                   ghi = ghi,
                                                   albedo = albedo,
                                                   surface_type = None,   # if provided, it overrides the albedo value
                                                   )

    # compute the direct (beam) component. Negative values (sun behind the plane) are zeroed
    beam_component = np.maximum(dni * np.cos(np.radians(aoi)), 0)

    # allocate the output variable and store every component
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


def inv_Pout_from_Pin(v_eff1,
                      v_eff2,
                      v_eff3,
                      p_in_v_eff1,
                      p_in_v_eff2,
                      p_in_v_eff3,
                      p_out_v_eff1,
                      p_out_v_eff2,
                      p_out_v_eff3,
                      p_dc_inv_in,
                      v_dc_in,
                      ):
    """
    Computes the inverter AC output power from its DC input power and MPPT input voltage, using the three
    efficiency curves stored in the OND file, expressed as Pout = f(Vin, Pin).

    The calculation is a double interpolation: first Pout is interpolated over Pin along each of the three
    voltage curves, then the three results are interpolated over the actual MPPT voltage.

    Args:
        v_eff1 (float):              [V] inverter input voltage at which p_in_v_eff1 / p_out_v_eff1 are known
        v_eff2 (float):              [V] inverter input voltage at which p_in_v_eff2 / p_out_v_eff2 are known
        v_eff3 (float):              [V] inverter input voltage at which p_in_v_eff3 / p_out_v_eff3 are known
        p_in_v_eff1 (np.ndarray):    [W] inverter DC input power points of the efficiency curve at v_eff1
        p_in_v_eff2 (np.ndarray):    [W] inverter DC input power points of the efficiency curve at v_eff2
        p_in_v_eff3 (np.ndarray):    [W] inverter DC input power points of the efficiency curve at v_eff3
        p_out_v_eff1 (np.ndarray):   [W] inverter AC output power resulting from p_in_v_eff1 at v_eff1
        p_out_v_eff2 (np.ndarray):   [W] inverter AC output power resulting from p_in_v_eff2 at v_eff2
        p_out_v_eff3 (np.ndarray):   [W] inverter AC output power resulting from p_in_v_eff3 at v_eff3
        p_dc_inv_in (pd.Series):     [W] actual DC input power of the inverter, for each simulation timestep
        v_dc_in (pd.Series):         [V] actual MPPT input voltage of the inverter, for each simulation timestep

    Returns:
        p_out_inv (pd.Series):  [W] AC output power of the inverter, for each simulation timestep
    """

    # create one Pout = f(Pin) interpolation function per voltage level
    interp_veff1 = interp1d(p_in_v_eff1, p_out_v_eff1, kind='quadratic', fill_value='extrapolate')
    interp_veff2 = interp1d(p_in_v_eff2, p_out_v_eff2, kind='quadratic', fill_value='extrapolate')
    interp_veff3 = interp1d(p_in_v_eff3, p_out_v_eff3, kind='quadratic', fill_value='extrapolate')

    # apply the interpolation functions above to the actual input power
    p_out1 = interp_veff1(p_dc_inv_in)
    p_out2 = interp_veff2(p_dc_inv_in)
    p_out3 = interp_veff3(p_dc_inv_in)

    # stack the three results so that they can be interpolated across voltages
    v_eff_values = [v_eff1, v_eff2, v_eff3]
    p_out_values = np.vstack([p_out1, p_out2, p_out3]).T

    # preallocate the output array
    p_out_final = np.zeros(len(v_dc_in))

    # interpolate over the voltage, timestep by timestep. Outside the known voltage range, the closest curve is used
    for i, v_eff in enumerate(v_dc_in):
        if v_eff <= v_eff1:
            p_out_final[i] = p_out_values[i, 0]
        elif v_eff >= v_eff3:
            p_out_final[i] = p_out_values[i, 2]
        else:
            interp_voltages = interp1d(v_eff_values, p_out_values[i, :], kind='quadratic')
            p_out_final[i] = interp_voltages(v_eff)

    # return the result with the same index as the input Series
    return pd.Series(p_out_final, index=p_dc_inv_in.index)


def ac_ohm_loss(p_ac,
                v_p2p,
                leng,
                material,
                section,
                ):
    """
    Computes the ohmic loss of a three phase AC cable, as a fraction of the power flowing through it.

    Args:
        p_ac (pd.Series or float):  [VA] apparent power flowing through ONE cable
        v_p2p (float):              [V] phase-to-phase voltage
        leng (float):               [m] cable length
        material (str):             [-] cable material, either 'cu'/'copper' or 'alu'/'aluminum'. Not case sensitive
        section (float):            [mm2] cable cross section

    Returns:
        ac_loss (pd.Series or float):  [-] ohmic loss as a unitless fraction in the range [0, 1]
    """

    # compute the phase current
    i_p2g = p_ac/math.sqrt(3)/v_p2p

    # resistivity as defined in PVsyst, in ohm*m
    if str.lower(str(material)) in ['cu', 'copper']:
        resistivity = Constants.resistivity_cu
    elif str.lower(str(material)) in ['alu', 'aluminum']:
        resistivity = Constants.resistivity_alu
    else:
        raise ValueError(f'Invalid cable material "{material}". Valid materials are copper, aluminum')

    # convert the section from mm2 to m2, then compute the cable resistance
    resistance = resistivity / (section/np.power(1000, 2)) * leng

    # three phase loss, expressed as a fraction of the power flowing through the cable
    ac_loss = 3 * resistance * np.power(i_p2g, 2) / p_ac

    return ac_loss


def vf_row_ground_2d_weights(standard_angles,
                             gcr,
                             pvrow_height,
                             pitch,
                             npoints=None,
                             n_pvrow=None,
                             max_rows=None,
                             ):
    """
    Computes how much each point of the ground contributes to the view of the ground from a row face.

    PVlib only exposes vf_row_ground_2d_integ, which is the view factor from the row to the WHOLE visible
    ground, already integrated. That is enough to know how much ground-reflected light reaches the modules
    when the ground is uniformly lit, but not when it is not: the strip right under the rows is permanently
    shaded while the strip halfway between rows is bright, and a row does not look at the two equally.
    This function returns the missing piece, i.e. the SHAPE of that view, so that the brightness of the
    ground can be averaged with the proper weights instead of uniformly over the pitch.

    The face is discretised into n_pvrow points. From each of them the contribution of a ground element is
    the 2D differential view factor between a surface element and a strip

        dF = ( sin(theta_2) - sin(theta_1) ) / 2

    with the angles measured from the normal of the face. A ground element is discarded when the line of
    sight from the row point to it is intercepted by one of the neighbouring rows.

    The geometry follows the same convention as PVlib's vf_ground_sky_2d: rows infinitely long, regularly
    spaced by 'pitch', with the row rotated so that its RIGHT edge is the higher one. The result only
    depends on the geometry, so it is computed on a few structure tilt angles and interpolated afterwards.

    Front and rear are returned together because they share the whole geometry and differ only by the
    direction of the normal, so computing them in one pass costs almost nothing more than one of the two.

    NOTE: the returned weights are normalised to sum to one along the ground, i.e. they carry the SHAPE of
    the view and not its magnitude. The magnitude keeps coming from PVlib's vf_row_ground_2d_integ, so that
    any first order error in the absolute value of these weights cancels out.

    Args:
        standard_angles (np.ndarray):  [deg] tilt angles of the STRUCTURE on which the geometry is computed,
                                       from 0 (flat) to 90. A tracker needs the whole range it sweeps,
                                       a fixed structure a single value
        gcr (float):                   [-] ground coverage ratio, row width divided by pitch
        pvrow_height (float):          [m] height of the centre of the PV row above the ground
        pitch (float):                 [m] distance between the axis of two adjacent rows
        npoints (int):                 [-] number of ground elements across one pitch. If None, Constants is used
        n_pvrow (int):                 [-] number of points the row face is discretised into. If None, Constants is used
        max_rows (int):                [-] neighbouring rows considered on each side. If None, Constants is used

    Returns:
        x_ground (np.ndarray):      [-] centre of each ground element, as a fraction of the pitch, shape (npoints,).
                                    x = 0 is the point right below the centre of the row, positive towards the right
        weights_front (np.ndarray): [-] contribution of each ground element to the view of the FRONT face,
                                    shape (len(standard_angles), npoints), each row summing to one
        weights_back (np.ndarray):  [-] same for the REAR face, shape (len(standard_angles), npoints)
    """

    if npoints is None:
        npoints = Constants.ground_npoints
    if n_pvrow is None:
        n_pvrow = Constants.ground_n_pvrow
    if max_rows is None:
        max_rows = Constants.ground_max_rows

    angles = np.radians(np.atleast_1d(np.asarray(standard_angles, dtype=float)))
    half_width = gcr*pitch/2

    # ---- the rows, as segments in the cross-axis plane, one per angle ----
    # the right edge is the higher one, as in PVlib's vf_ground_sky_2d
    edge_right = np.stack([+half_width*np.cos(angles), pvrow_height + half_width*np.sin(angles)], axis=-1)
    edge_left = np.stack([-half_width*np.cos(angles), pvrow_height - half_width*np.sin(angles)], axis=-1)

    # points along the face, avoiding the very edges where the view factor is singular.  shape (angles, n_pvrow, 2)
    fraction = np.linspace(0.05, 0.95, n_pvrow)
    row_points = edge_left[:, np.newaxis, :] + fraction[np.newaxis, :, np.newaxis]*(edge_right - edge_left)[:, np.newaxis, :]

    # normal of each face. With the right edge higher, the front looks to the left and the rear to the right
    normal_front = np.stack([-np.sin(angles), np.cos(angles)], axis=-1)
    normals = np.stack([normal_front, -normal_front], axis=0)                  # shape (2, angles, 2), front then rear
    tangents = np.stack([-normals[..., 1], normals[..., 0]], axis=-1)

    # ---- the ground, discretised across one pitch and centred below the row ----
    edges = np.linspace(-0.5, 0.5, npoints+1)*pitch
    x_centre = (edges[:-1] + edges[1:])/2

    # ---- view factor from every point of every face to every ground element ----
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

    # ---- lines of sight intercepted by the neighbouring rows ----
    # the blockers are the same rows shifted by a whole number of pitches
    shift = np.array([k*pitch for k in range(-max_rows, max_rows+1) if k != 0])
    blocker_1 = edge_left[:, np.newaxis, :] + np.stack([shift, np.zeros_like(shift)], axis=-1)[np.newaxis, :, :]
    blocker_2 = edge_right[:, np.newaxis, :] + np.stack([shift, np.zeros_like(shift)], axis=-1)[np.newaxis, :, :]

    # segment from a row point P to a ground centre C, against the blocker segment Q1-Q2.
    # Shapes are broadcast to (angles, n_pvrow, npoints, blockers)
    p_x = row_points[..., 0][:, :, np.newaxis, np.newaxis]
    p_y = row_points[..., 1][:, :, np.newaxis, np.newaxis]
    c_x = x_centre[np.newaxis, np.newaxis, :, np.newaxis]
    c_y = 0.0
    q1_x = blocker_1[..., 0][:, np.newaxis, np.newaxis, :]
    q1_y = blocker_1[..., 1][:, np.newaxis, np.newaxis, :]
    q2_x = blocker_2[..., 0][:, np.newaxis, np.newaxis, :]
    q2_y = blocker_2[..., 1][:, np.newaxis, np.newaxis, :]

    def _side(origin_x, origin_y, a_x, a_y, b_x, b_y):
        """Cross product telling on which side of the segment origin-a the point b lies."""
        return (a_x-origin_x)*(b_y-origin_y) - (a_y-origin_y)*(b_x-origin_x)

    side_1 = _side(p_x, p_y, c_x, c_y, q1_x, q1_y)
    side_2 = _side(p_x, p_y, c_x, c_y, q2_x, q2_y)
    side_3 = _side(q1_x, q1_y, q2_x, q2_y, p_x, p_y)
    side_4 = _side(q1_x, q1_y, q2_x, q2_y, c_x, c_y)

    blocked = np.any((side_1*side_2 < 0) & (side_3*side_4 < 0), axis=-1)   # shape (angles, n_pvrow, npoints)

    # ---- sum over the points of the face, then normalise along the ground ----
    weights = np.sum(np.where(visible & ~blocked[np.newaxis, ...], delta_vf, 0.0), axis=2)   # shape (2, angles, npoints)

    total = np.sum(weights, axis=-1, keepdims=True)
    weights = np.divide(weights, total, out=np.zeros_like(weights), where=total > 0)

    return x_centre/pitch, weights[0], weights[1]


def ground_brightness_2d(x,
                         surface_tilt,
                         gcr,
                         pvrow_height,
                         pitch,
                         diffuse_fraction,
                         tan_psi,
                         max_rows=None,
                         ):
    """
    Computes the brightness of each point of the ground, relative to an open, unobstructed ground.

    The ground between the rows is not uniformly lit: it only sees part of the sky, and part of it is in the
    shadow the rows cast. This function returns both effects combined, point by point:

        brightness(x) = diffuse_fraction * (sky seen from x) + (1 - diffuse_fraction) * (x is lit by the sun)

    Args:
        x (np.ndarray):             [-] ground positions, as a fraction of the pitch, with the same
                                    convention as vf_row_ground_2d_weights
        surface_tilt (float):       [deg] tilt of the face being evaluated, from 0 to 180. Only the tilt of
                                    the structure is taken from it, as the ground does not depend on the face
        gcr (float):                [-] ground coverage ratio
        pvrow_height (float):       [m] height of the centre of the PV row above the ground
        pitch (float):              [m] distance between the axis of two adjacent rows
        diffuse_fraction (float):   [-] share of the horizontal irradiance that is diffuse
        tan_psi (float):            [-] tangent of the solar projection angle on the cross-axis plane,
                                    positive when the sun is on the side the FRONT of the modules faces
        max_rows (int):             [-] neighbouring rows considered. If None, Constants is used

    Returns:
        np.ndarray:  [-] brightness of each ground point, between 0 and 1
    """

    if max_rows is None:
        max_rows = Constants.ground_max_rows_sky

    pvrow_tilt_deg = min(surface_tilt, 180-surface_tilt)
    pvrow_tilt = np.radians(pvrow_tilt_deg)
    half_width = gcr*pitch/2

    # DIFFUSE PART - how much sky each ground point still sees
    vf_sky = utils.vf_ground_sky_2d(rotation = pvrow_tilt_deg,
                                    gcr = gcr,
                                    x = x,
                                    pitch = pitch,
                                    height = pvrow_height,
                                    max_rows = max_rows,
                                    )
    vf_sky = np.asarray(vf_sky).ravel()

    # BEAM PART - where the shadow of the rows falls.
    # With the right edge of the row higher, the front of the modules looks to the left, so a sun on the
    # front side casts the shadow towards the right: the horizontal shift per unit of height is +tan_psi
    edge_right = np.array([+half_width*np.cos(pvrow_tilt), pvrow_height + half_width*np.sin(pvrow_tilt)])
    edge_left = np.array([-half_width*np.cos(pvrow_tilt), pvrow_height - half_width*np.sin(pvrow_tilt)])

    shadow_1 = (edge_right[0] + edge_right[1]*tan_psi)/pitch
    shadow_2 = (edge_left[0] + edge_left[1]*tan_psi)/pitch

    shadow_centre = (shadow_1 + shadow_2)/2
    shadow_half_width = abs(shadow_1 - shadow_2)/2

    # the shadow pattern repeats every pitch, hence the distance is taken modulo one
    distance_to_shadow = np.abs(((x - shadow_centre + 0.5) % 1.0) - 0.5)
    lit = (distance_to_shadow > shadow_half_width).astype(float)

    return diffuse_fraction*vf_sky + (1-diffuse_fraction)*lit
