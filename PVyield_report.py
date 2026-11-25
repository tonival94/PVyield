# -*- coding: utf-8 -*-
"""
Support functions for PVyield: export of the simulation inputs and results.

The main entry point is post_processing(), which writes an Excel workbook holding:

    SHEET 1 - 'Simulation'
        Every input parameter of the simulation, with its unit of measurement and a short note, split into
        a block common to the whole project and a block holding one column per PV field.
        The loss chain and the Performance Ratio of the PROJECT follow at the bottom of the sheet.

    SHEET 2 - 'Hourly energy'
        The net energy delivered at the Point of Connection for each of the 8760 hours of the year,
        plus a monthly summary computed with live Excel formulas.

@author: Antonio Vallone (-T O N I-)

Libraries requirements:
> openpyxl==3.1.5
"""

from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

from PVyield_classes import Results, WarningLog

_CURRENT_DIR = Path(__file__).parent.resolve()

# ---------------------------------------------------------------------------
#   FORMATTING CONSTANTS
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
#   PARAMETER TABLES
# ---------------------------------------------------------------------------
# Each row is a tuple: (section, parameter name, value, unit, note).
# The note explains what the parameter drives in the simulation, so that the workbook can be read
# by someone who has never seen the launcher.


def _format_value(value):
    """
    Converts a parameter value into something Excel can store.

    Args:
        value (any):  the raw value read from the input classes

    Returns:
        str, float, int or None:  a value openpyxl is able to write
    """

    if value is None:
        return 'not used'

    if isinstance(value, bool):
        return 'Yes' if value else 'No'

    if isinstance(value, (tuple, list, np.ndarray)):
        # short arrays are written inline, long ones are only summarised
        values = list(np.ravel(value))
        if len(values) > 12:
            return f'{len(values)} points (see source file)'
        return ', '.join(f'{v:g}' if isinstance(v, (int, float, np.number)) else str(v) for v in values)

    if isinstance(value, (np.integer, np.floating)):
        return value.item()

    return value


def _project_rows(pv_project):
    """
    Builds the table of the input parameters common to the whole project.

    Args:
        pv_project (PV_Project):  PV_Project class object

    Returns:
        list[tuple]:  rows shaped as (section, parameter, value, unit, note)
    """

    rows = []

    rows.append(('Location', 'Latitude', pv_project.coordinates.latitude, 'deg', 'Decimal degrees, positive in the North hemisphere'))
    rows.append(('Location', 'Longitude', pv_project.coordinates.longitude, 'deg', 'Decimal degrees, positive in the East hemisphere'))

    rows.append(('TMY', 'TMY source', pv_project.tmy.TMY_source, '-', 'Source of the meteorological data, common to every PV field'))
    rows.append(('TMY', 'TMY file name', pv_project.tmy.TMY_file_name or 'not used', '-', 'Only used when the TMY is read from a file'))
    # NOTE: the API key is deliberately NOT exported, so that the workbook can be shared safely

    rows.append(('HV trafo', 'Rated power', pv_project.hv_trafo.p_nom, 'MVA', 'MV/HV transformer. If not used, the PoC sits on the MV network'))
    rows.append(('HV trafo', 'Iron loss', pv_project.hv_trafo.iron_loss, '%', 'No load loss of the MV/HV transformer'))
    rows.append(('HV trafo', 'Copper loss', pv_project.hv_trafo.copper_loss, '%', 'Load loss of the MV/HV transformer at rated power'))

    rows.append(('HV connection', 'Line length', pv_project.hv_connection.line_len, 'km', 'HV line from the MV/HV transformer to the Point of Connection'))

    rows.append(('Project curtailment', 'Application point', pv_project.project_curtailment.grid_curtail_point, '-', 'Node of the common part where the project curtailment is applied'))
    rows.append(('Project curtailment', 'Curtailment power', pv_project.project_curtailment.grid_curtail_power, 'MVA', 'AC power cap applied at the point above'))

    rows.append(('Availability', 'Grid availability', pv_project.grid_availability, '%', 'Share of the year the grid is able to accept the production'))

    rows.append(('Summary', 'Number of PV fields', len(pv_project.pv_fields), '-', 'PV fields sharing the same location and the same TMY'))
    rows.append(('Summary', 'Total DC power', sum(f.power.p_dc_peak for f in pv_project.pv_fields), 'MWp', 'Sum of the peak DC power of every PV field'))
    rows.append(('Summary', 'Total inverter power', sum(f.power.p_ac_nom for f in pv_project.pv_fields), 'MVA', 'Sum of the nominal inverter power of every PV field'))
    rows.append(('Summary', 'Field feeding the PoC', pv_project.exit_field.name, '-', 'PV field owning the last MV line of the tree'))

    return rows


def _field_rows(field):
    """
    Builds the table of the input parameters of one PV field.

    Every PV field produces the same rows in the same order, so that the tables of different fields can be
    put side by side as columns of a single table.

    Args:
        field (PV_Field):  PV_Field class object

    Returns:
        list[tuple]:  rows shaped as (section, parameter, value, unit, note)
    """

    rows = []

    # ----- POWER -----
    rows.append(('Power', 'Peak DC power', field.power.p_dc_peak, 'MWp', 'Rounded to an integer number of strings during the input check'))
    rows.append(('Power', 'Nominal inverter power', field.power.p_ac_nom, 'MVA', 'Rounded to an integer number of inverters during the input check'))
    rows.append(('Power', 'DC/AC ratio', field.power.p_dc_peak/field.power.p_ac_nom, '-', 'Peak DC power divided by the nominal inverter power'))

    # ----- MOUNTING STRUCTURE -----
    rows.append(('Mounting structure', 'Mount type', field.mms.mount_type, '-', 'Single axis tracker or fixed tilt'))
    rows.append(('Mounting structure', 'Module configuration', field.mms.mod_configuration, '-', 'Modules per table and orientation, P = portrait, L = landscape'))
    rows.append(('Mounting structure', 'Modules per string', field.mms.modules_per_string, '-', 'Drives the string voltage and the MPPT window'))
    rows.append(('Mounting structure', 'Pitch', field.mms.pitch, 'm', 'Row to row distance, drives the ground coverage ratio'))
    rows.append(('Mounting structure', 'Tilt angle', field.mms.tilt_angle, 'deg', 'Max rotation range for a tracker, surface tilt for a fixed structure'))
    rows.append(('Mounting structure', 'Axis azimuth', field.mms.axis_azimuth, 'deg', 'Clockwise from North: 90 = East, 180 = South, 270 = West'))
    rows.append(('Mounting structure', 'Ground gap', field.mms.ground_gap, 'm', 'Distance between ground and the lower edge of the table'))
    rows.append(('Mounting structure', 'Backtracking', field.mms.backtrack, '-', 'Tracker only. Avoids row to row shading at low sun elevation'))

    # ----- GROUND -----
    rows.append(('Ground', 'Albedo', field.ground.albedo, '-', 'Ground reflectivity. Overridden by the TMY when it carries albedo data'))
    rows.append(('Ground', 'Slope tilt', field.ground.slope_tilt, '%', 'Ground slope, drives the axis tilt and the cross axis tilt'))
    rows.append(('Ground', 'Slope azimuth', field.ground.slope_azimuth, 'deg', 'Direction where a ball would roll if left on the slope'))

    # ----- MODULE -----
    rows.append(('PV module', 'Manufacturer', field.module.manufacturer, '-', 'From the PAN file when PAN_import is active'))
    rows.append(('PV module', 'Model', field.module.model, '-', 'From the PAN file when PAN_import is active'))
    rows.append(('PV module', 'PAN file', field.module.PAN_file_name, '-', 'Source of the module parameters'))
    rows.append(('PV module', 'Nominal power', field.module.Pnom, 'Wp', 'Module DC power at STC'))
    rows.append(('PV module', 'Height', field.module.height, 'm', 'Drives the table width when mounted in portrait'))
    rows.append(('PV module', 'Width', field.module.width, 'm', 'Drives the table width when mounted in landscape'))
    rows.append(('PV module', 'Cells in series', field.module.cells_in_series, '-', 'Full cut equivalent number of cells'))
    rows.append(('PV module', 'Cells in parallel', field.module.cells_in_parallel, '-', '2 for half cut modules, 1 for full cut modules'))
    rows.append(('PV module', 'Bypass diodes', field.module.n_diodes, '-', 'Number of bypass diodes'))
    rows.append(('PV module', 'Short circuit current', field.module.I_L_ref, 'A', 'Photocurrent at STC, single diode model input'))
    rows.append(('PV module', 'Open circuit voltage', field.module.Voc, 'V', 'Drives the max number of modules per string'))
    rows.append(('PV module', 'MPP current', field.module.imp_ref, 'A', 'Current at maximum power point at STC'))
    rows.append(('PV module', 'MPP voltage', field.module.vmp_ref, 'V', 'Voltage at maximum power point at STC'))
    rows.append(('PV module', 'Isc temp coefficient', field.module.alpha_sc, 'mA/degC', 'Short circuit current temperature coefficient'))
    rows.append(('PV module', 'Diode ideality factor', field.module.gamma_ref, '-', 'PVsyst single diode model parameter'))
    rows.append(('PV module', 'Ideality temp coefficient', field.module.mu_gamma, '1/K', 'PVsyst single diode model parameter'))
    rows.append(('PV module', 'Shunt resistance', field.module.R_sh_ref, 'ohm', 'Shunt resistance at STC'))
    rows.append(('PV module', 'Shunt resistance at 0 irr', field.module.R_sh_0, 'ohm', 'Shunt resistance at null irradiance'))
    rows.append(('PV module', 'Shunt exponent', field.module.R_sh_exp, '-', 'PVsyst single diode model parameter'))
    rows.append(('PV module', 'Series resistance', field.module.R_s, 'ohm', 'Series resistance at STC'))
    rows.append(('PV module', 'Absorption coefficient', field.module.alpha_absorption, '-', 'Drives the cell temperature through the PVsyst thermal model'))
    rows.append(('PV module', 'Bifaciality factor', field.module.bifaciality, '-', 'Rear side efficiency relative to the front side'))
    rows.append(('PV module', 'Bypass diode voltage', field.module.V_bypass_diode, 'V', 'Voltage drop across a conducting bypass diode'))
    rows.append(('PV module', 'Substring layout', field.module.substring_type, '-', 'Half cut or full cut, drives the partial shading loss'))

    # ----- INVERTER -----
    rows.append(('Inverter', 'Manufacturer', field.inverter.manufacturer, '-', 'From the OND file when OND_import is active'))
    rows.append(('Inverter', 'Model', field.inverter.model, '-', 'From the OND file when OND_import is active'))
    rows.append(('Inverter', 'OND file', field.inverter.OND_file_name, '-', 'Source of the inverter parameters'))
    rows.append(('Inverter', 'Nominal AC power', field.inverter.Pnom_ac, 'kW', 'Continuous AC output power of a single inverter'))
    rows.append(('Inverter', 'Max AC power', field.inverter.P_max_ac, 'kW', 'Peak AC output power of a single inverter'))
    rows.append(('Inverter', 'Number of inverters', int(np.ceil(field.power.p_ac_nom*1000/field.inverter.Pnom_ac)), '-', 'Derived from the nominal inverter power of the field'))
    rows.append(('Inverter', 'Min MPPT voltage', field.inverter.v_min_mppt, 'V', 'Lower boundary of the MPPT voltage scan'))
    rows.append(('Inverter', 'Max MPPT voltage', field.inverter.v_max_mppt, 'V', 'Upper boundary of the MPPT voltage scan'))
    rows.append(('Inverter', 'Max MPPT current', field.inverter.i_max_mppt, 'A', 'Current limit of each MPPT input'))
    rows.append(('Inverter', 'String inputs', field.inverter.n_string_inputs, '-', 'Max number of strings connected to one inverter'))
    rows.append(('Inverter', 'Number of MPPT', field.inverter.n_mppt, '-', 'Drives how the strings are distributed inside the inverter'))
    rows.append(('Inverter', 'AC output voltage', field.inverter.v_out, 'V', 'Drives the LV cable ohmic loss'))
    rows.append(('Inverter', 'Temp at nominal power', field.inverter.T_p_nom, 'degC', 'Above this temperature the inverter derates'))
    rows.append(('Inverter', 'Temp at max power', field.inverter.T_p_max, 'degC', 'Up to this temperature the inverter delivers its max power'))

    # ----- LV/MV TRAFO -----
    rows.append(('LV/MV trafo', 'Output voltage', field.mv_trafo.v_out, 'kV', 'MV network voltage, common to every PV field'))
    rows.append(('LV/MV trafo', 'Max rated power', field.mv_trafo.p_max, 'MVA', 'Drives the number of transformers in the field'))
    rows.append(('LV/MV trafo', 'Iron loss', field.mv_trafo.iron_loss, '%', 'No load loss'))
    rows.append(('LV/MV trafo', 'Copper loss', field.mv_trafo.copper_loss, '%', 'Load loss at rated power'))

    # ----- LOSSES -----
    rows.append(('Losses', 'IAM model', field.loss.iam_model, '-', 'Incidence angle modifier model'))
    rows.append(('Losses', 'Soiling loss', field.loss.soiling_loss, '%', 'Dust and dirt on the module front glass'))
    rows.append(('Losses', 'Rear shade factor', field.loss.rear_shade_factor, '%', 'Light blocked by the structure behind the modules'))
    rows.append(('Losses', 'Thermal loss factor', field.loss.thermal_loss_factor, 'W/m2/K', 'PVsyst thermal model coefficient'))
    rows.append(('Losses', 'Module quality', field.loss.module_quality, '%', 'A negative value is a gain'))
    rows.append(('Losses', 'LID + LeTID loss', field.loss.lid_letid_loss, '%', 'Light and elevated temperature induced degradation'))
    rows.append(('Losses', 'First year degradation', field.loss.first_year_degrad, '%', 'Applied before the curtailment'))
    rows.append(('Losses', 'Mismatch front', field.loss.mismatch_loss_front, '%', 'Array mismatch on the front side'))
    rows.append(('Losses', 'Mismatch back', field.loss.mismatch_loss_back, '%', 'Array mismatch on the rear side'))
    rows.append(('Losses', 'DC ohmic loss at STC', field.loss.dc_ohmic_loss_STC, '%', 'String cabling loss at STC, scaled with the actual current'))
    rows.append(('Losses', 'Auxiliaries loss', field.loss.auxiliaries_loss, '%', 'Trackers, SCADA, lighting, HVAC of the cabins'))
    rows.append(('Losses', 'PV availability', field.loss.pv_availability, '%', 'Share of the year the field is able to produce'))
    rows.append(('Losses', 'LV cable length', field.loss.ac_lv_cable_len, 'm', 'Average length from inverter to transformer'))
    rows.append(('Losses', 'LV cable section', field.loss.ac_lv_cable_sect, 'mm2', 'Average section of the LV cable'))
    rows.append(('Losses', 'LV cable material', field.loss.ac_lv_cable_material, '-', 'Copper or aluminium'))
    rows.append(('Losses', 'MV cable length', field.loss.ac_mv_cable_len, 'm', 'Average length from transformer to the field MV cabin'))
    rows.append(('Losses', 'MV cable section', field.loss.ac_mv_cable_sect, 'mm2', 'Average section of the MV cable'))
    rows.append(('Losses', 'MV cable material', field.loss.ac_mv_cable_material, '-', 'Copper or aluminium'))

    # ----- FIELD CURTAILMENT -----
    rows.append(('Field curtailment', 'Application point', field.field_curtailment.grid_curtail_point, '-', 'Inverter output or field MV cabin'))
    rows.append(('Field curtailment', 'Curtailment power', field.field_curtailment.grid_curtail_power, 'MVA', 'AC power cap of the whole field at the point above'))

    # ----- MV CONNECTION -----
    rows.append(('MV connection', 'Destination', field.mv_connection.destination_field, '-', 'Node the MV line of this field is connected to'))
    rows.append(('MV connection', 'Line length', field.mv_connection.line_len, 'km', 'Length of the MV line leaving this field'))
    rows.append(('MV connection', 'Number of cables', field.mv_connection.line_nr, '-', 'Parallel cables of the MV line'))
    rows.append(('MV connection', 'Cable section', field.mv_connection.line_sect, 'mm2', 'Section of each cable of the MV line'))
    rows.append(('MV connection', 'Cable material', field.mv_connection.material, '-', 'Copper or aluminium'))
    rows.append(('MV connection', 'Automatic sizing', field.mv_connection.auto_sizing, '-', 'Not implemented yet, the size must be declared manually'))

    return rows


# ---------------------------------------------------------------------------
#   WRITING HELPERS
# ---------------------------------------------------------------------------

def _write_title(sheet, row, text, n_columns):
    """Writes a full width title bar and returns the next free row."""

    sheet.cell(row=row, column=1, value=text)
    for col in range(1, n_columns+1):
        sheet.cell(row=row, column=col).fill = _TITLE_FILL
        sheet.cell(row=row, column=col).font = _TITLE_FONT
    sheet.row_dimensions[row].height = 22
    return row + 1


def _write_section(sheet, row, text, n_columns):
    """Writes a section bar and returns the next free row."""

    sheet.cell(row=row, column=1, value=text)
    for col in range(1, n_columns+1):
        sheet.cell(row=row, column=col).fill = _SECTION_FILL
        sheet.cell(row=row, column=col).font = _SECTION_FONT
    return row + 1


def _write_header(sheet, row, labels):
    """Writes a table header and returns the next free row."""

    for col, label in enumerate(labels, start=1):
        cell = sheet.cell(row=row, column=col, value=label)
        cell.font = _HEADER_FONT
        cell.fill = _HEADER_FILL
        cell.border = _THIN_BORDER
        cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
    return row + 1


# ---------------------------------------------------------------------------
#   SHEET BUILDERS
# ---------------------------------------------------------------------------

def _build_simulation_sheet(sheet, pv_project, RESULTS):
    """
    Fills the first sheet with the simulation inputs and the project results.

    Args:
        sheet (openpyxl worksheet):  destination sheet
        pv_project (PV_Project):     PV_Project class object
        RESULTS (pd.DataFrame):      RESULTS DataFrame returned by PVyield
    """

    field_names = [f.name for f in pv_project.pv_fields]
    # The results block is the widest one: Loss name, Unit, one column per PV field, plus the project column.
    # The notes of the input blocks sit one column further right, so that they never share a column
    # with the project results and each column can be given its own width.
    note_column = 4 + len(field_names)
    n_columns = note_column

    row = _write_title(sheet, 1, 'PVyield  -  simulation inputs and results', n_columns)
    sheet.cell(row=row, column=1, value=f'Generated on {datetime.now():%Y-%m-%d %H:%M}').font = _NOTE_FONT
    row += 2

    # =========================================================================
    # PROJECT INPUTS
    # =========================================================================
    row = _write_section(sheet, row, 'PROJECT INPUTS  -  common to every PV field', n_columns)
    row = _write_header(sheet, row, ['Parameter', 'Unit', 'Value'] + ['']*(note_column-4) + ['Note'])   # blanks pad the header up to the note column

    current_section = None
    for section, parameter, value, unit, note in _project_rows(pv_project):

        # a blank line separates the groups of parameters
        if current_section is not None and section != current_section:
            row += 1
        current_section = section

        sheet.cell(row=row, column=1, value=f'{section}  |  {parameter}').font = _BODY_FONT
        sheet.cell(row=row, column=2, value=unit).font = _BODY_FONT
        sheet.cell(row=row, column=2).alignment = Alignment(horizontal='center')

        cell = sheet.cell(row=row, column=3, value=_format_value(value))
        cell.font = _BODY_FONT
        if isinstance(cell.value, float):
            cell.number_format = _NUM_FORMAT

        sheet.cell(row=row, column=note_column, value=note).font = _NOTE_FONT
        row += 1

    row += 2

    # =========================================================================
    # PV FIELD INPUTS
    # =========================================================================
    row = _write_section(sheet, row, 'PV FIELD INPUTS  -  one column per PV field', n_columns)
    row = _write_header(sheet, row, ['Parameter', 'Unit'] + field_names + [''] + ['Note'])   # the blank keeps the note aligned on note_column

    # every field produces the same rows in the same order, so they can be zipped together
    field_tables = [_field_rows(f) for f in pv_project.pv_fields]

    current_section = None
    for row_index, reference_row in enumerate(field_tables[0]):
        section, parameter, _, unit, note = reference_row

        if current_section is not None and section != current_section:
            row += 1
        current_section = section

        sheet.cell(row=row, column=1, value=f'{section}  |  {parameter}').font = _BODY_FONT
        sheet.cell(row=row, column=2, value=unit).font = _BODY_FONT
        sheet.cell(row=row, column=2).alignment = Alignment(horizontal='center')

        for field_index, table in enumerate(field_tables):
            cell = sheet.cell(row=row, column=3+field_index, value=_format_value(table[row_index][2]))
            cell.font = _BODY_FONT
            cell.alignment = Alignment(horizontal='center')
            if isinstance(cell.value, float):
                cell.number_format = _NUM_FORMAT

        sheet.cell(row=row, column=note_column, value=note).font = _NOTE_FONT
        row += 1

    row += 2

    # =========================================================================
    # RESULTS
    # =========================================================================
    # Differently from the console print, which only shows the project column to keep the table readable,
    # the Excel report holds the full loss chain: one column per PV field, plus the project column.
    row = _write_section(sheet, row, 'RESULTS  -  loss chain of each PV field and of the whole project', n_columns)
    row = _write_header(sheet, row, ['Loss / gain', 'Unit'] + field_names + [Results.PROJECT_COLUMN])

    highlight_rows = [Results.PR_PV_ONLY, Results.YIELD_PV_ONLY,
                      Results.PR_PV_MV, Results.YIELD_PV_MV,
                      Results.PR_PROJECT, Results.YIELD_PROJECT]

    # the value columns of the RESULTS DataFrame, in the same order as the header written above
    result_columns = field_names + [Results.PROJECT_COLUMN]
    last_result_column = 2 + len(result_columns)

    for label in RESULTS.index:

        # the spacers of the loss chain become blank separator lines
        if label in Results.SPACERS:
            row += 1
            continue

        is_highlight = label in highlight_rows
        font = _BOLD_FONT if is_highlight else _BODY_FONT

        sheet.cell(row=row, column=1, value=label).font = font
        sheet.cell(row=row, column=2, value=RESULTS.loc[label, 'Units']).font = font
        sheet.cell(row=row, column=2).alignment = Alignment(horizontal='center')

        for column_index, result_column in enumerate(result_columns):
            value = RESULTS.loc[label, result_column]
            # rows that are meaningless for a single PV field hold pd.NA and are left blank
            cell = sheet.cell(row=row, column=3+column_index, value=None if pd.isna(value) else float(value))
            cell.font = font
            cell.number_format = _NUM_FORMAT
            cell.alignment = Alignment(horizontal='center')

        if is_highlight:
            for col in range(1, last_result_column+1):
                sheet.cell(row=row, column=col).fill = _HIGHLIGHT_FILL

        row += 1

    row += 1
    sheet.cell(row=row, column=1, value='Note: a negative value is a loss, a positive value is a gain. '
                                        'Each Performance Ratio is the running product of (1 + value) over the rows above it. '
                                        'Blank cells in the PV field columns hold values that only make sense at project level '
                                        '(project curtailment, HV parts and grid availability).').font = _NOTE_FONT

    # column widths
    sheet.column_dimensions['A'].width = 46
    sheet.column_dimensions['B'].width = 10
    # one column per PV field, plus the project column of the results block
    for index in range(len(field_names)+1):
        sheet.column_dimensions[get_column_letter(3+index)].width = 22
    sheet.column_dimensions[get_column_letter(note_column)].width = 70

    sheet.freeze_panes = 'C4'


def _build_energy_sheet(sheet, energy_hourly):
    """
    Fills the second sheet with the hourly net energy and a monthly summary.

    Args:
        sheet (openpyxl worksheet):   destination sheet
        energy_hourly (pd.Series):    [kWh] hourly net energy at the Point of Connection
    """

    _write_title(sheet, 1, 'PVyield  -  hourly net energy at the Point of Connection', 6)

    header_row = 3
    _write_header(sheet, header_row, ['Timestamp (local time)', 'Month', 'Net energy [kWh]'])

    # Excel cannot store a timezone aware timestamp, so the timezone is removed.
    # The values are already expressed in the local time of the site
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

    # ---------------------------------------------------------------------
    # MONTHLY SUMMARY
    # ---------------------------------------------------------------------
    # Built with live formulas, so that the summary follows any edit made to the hourly values

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

        # SUMIFS is an Excel 2007 function, so it is evaluated correctly by every spreadsheet application
        formula = (f'=SUMIFS($C${first_data_row}:$C${last_data_row},'
                   f'$B${first_data_row}:$B${last_data_row},{month_index+1})/1000')
        cell = sheet.cell(row=current_row, column=summary_column+1, value=formula)
        cell.font = _BODY_FONT
        cell.number_format = _ENERGY_FORMAT

    total_row = first_data_row + 12
    cell = sheet.cell(row=total_row, column=summary_column, value='TOTAL')
    cell.font = _BOLD_FONT
    cell.fill = _HIGHLIGHT_FILL

    cell = sheet.cell(row=total_row, column=summary_column+1,
                      value=f'=SUM({get_column_letter(summary_column+1)}{first_data_row}:'
                            f'{get_column_letter(summary_column+1)}{total_row-1})')
    cell.font = _BOLD_FONT
    cell.fill = _HIGHLIGHT_FILL
    cell.number_format = _ENERGY_FORMAT

    sheet.column_dimensions['A'].width = 22
    sheet.column_dimensions['B'].width = 8
    sheet.column_dimensions['C'].width = 16
    sheet.column_dimensions['D'].width = 4
    sheet.column_dimensions['E'].width = 14
    sheet.column_dimensions['F'].width = 18

    sheet.freeze_panes = f'A{first_data_row}'


def _build_warnings_sheet(sheet, warning_messages):
    """
    Fills the third sheet with the warnings raised while the simulation was running, one per row.

    Warnings never stop the simulation: they report a missing or inconsistent input that has been replaced
    by a standard value, or a result worth a second look. Having them stored next to the results makes it
    possible to tell, months later, which numbers of the report rest on an assumption.

    Args:
        sheet (openpyxl worksheet):     destination sheet
        warning_messages (list[dict]):  warnings collected by the WarningLog class, each holding
                                        the keys 'category', 'message' and 'source'
    """

    _write_title(sheet, 1, 'PVyield  -  warnings raised during the simulation', 4)

    if not warning_messages:
        cell = sheet.cell(row=3, column=1, value='No warning was raised: every input was read as provided, '
                                                 'with no standard value assumed in its place.')
        cell.font = _BODY_FONT
        sheet.column_dimensions['A'].width = 100
        return

    sheet.cell(row=2, column=1, value=f'{len(warning_messages)} warning(s) raised. Each one reports an input that was '
                                      f'missing or inconsistent and has been replaced by a standard value, or a result '
                                      f'worth checking. The simulation completed regardless.').font = _NOTE_FONT

    header_row = 4
    _write_header(sheet, header_row, ['#', 'Type', 'Warning message', 'Raised in'])

    for index, warning in enumerate(warning_messages):
        current_row = header_row + 1 + index

        cell = sheet.cell(row=current_row, column=1, value=index+1)
        cell.font = _BODY_FONT
        cell.alignment = Alignment(horizontal='center', vertical='top')

        cell = sheet.cell(row=current_row, column=2, value=warning.get('category', ''))
        cell.font = _BODY_FONT
        cell.alignment = Alignment(horizontal='center', vertical='top')

        # the message is written on a single row, wrapped, so that one warning always equals one row
        cell = sheet.cell(row=current_row, column=3, value=warning.get('message', ''))
        cell.font = _BODY_FONT
        cell.alignment = Alignment(vertical='top', wrap_text=True)

        cell = sheet.cell(row=current_row, column=4, value=warning.get('source', ''))
        cell.font = _NOTE_FONT
        cell.alignment = Alignment(horizontal='left', vertical='top')

    sheet.column_dimensions['A'].width = 6
    sheet.column_dimensions['B'].width = 16
    sheet.column_dimensions['C'].width = 110
    sheet.column_dimensions['D'].width = 26

    sheet.freeze_panes = f'A{header_row+1}'


# ---------------------------------------------------------------------------
#   FILE NAME
# ---------------------------------------------------------------------------

def _build_file_name(pv_project):
    """
    Builds a self explanatory file name, and makes sure it does not overwrite an existing file.

    Args:
        pv_project (PV_Project):  PV_Project class object

    Returns:
        Path:  full path of the file to be written
    """

    p_dc = sum(f.power.p_dc_peak for f in pv_project.pv_fields)
    p_ac = sum(f.power.p_ac_nom for f in pv_project.pv_fields)

    name = (f'{datetime.now():%Y%m%d_%H%M}'
            f'_{p_dc:.2f}MWp'
            f'_{p_ac:.2f}MVA'
            f'_{len(pv_project.pv_fields)}fields')

    # remove the characters that are not allowed in a file name
    for character in '\\/:*?"<>|':
        name = name.replace(character, '-')

    results_dir = _CURRENT_DIR / 'Results'
    results_dir.mkdir(exist_ok=True)

    file_path = results_dir / f'{name}.xlsx'

    # if a file with the same name already exists, a progressive number is appended
    counter = 1
    while file_path.exists():
        file_path = results_dir / f'{name}_({counter}).xlsx'
        counter += 1

    return file_path


# ---------------------------------------------------------------------------
#   MAIN ENTRY POINT
# ---------------------------------------------------------------------------

def excel_generator(pv_project, RESULTS, energy_hourly, warning_messages=None, file_name=None):
    """
    Exports the simulation inputs and results into an Excel workbook.

    The workbook holds three sheets:
        'Simulation'    every input parameter with its unit and a short note, plus the loss chain of each
                        PV field and of the whole project
        'Hourly energy' the net energy delivered at the Point of Connection for each hour of the year,
                        plus a monthly summary computed with live formulas
        'Warnings'      every warning raised while the simulation was running, one per row

    Args:
        pv_project (PV_Project):            PV_Project class object used to run the simulation
        RESULTS (pd.DataFrame):             RESULTS DataFrame returned by PVyield
        energy_hourly (pd.Series):          [kWh] hourly net energy returned by PVyield
        warning_messages (list, optional):  warnings to be listed in the 'Warnings' sheet. If None, they are
                                            read from the WarningLog class, which collects them automatically
                                            as soon as WarningLog.start() has been called by the launcher
        file_name (str, optional):          name of the output file. If None, a self explanatory name is built
                                            automatically and the file is saved in the 'Results' subfolder

    Returns:
        Path:  full path of the workbook that has been written
    """

    if file_name is None:
        file_path = _build_file_name(pv_project)
    else:
        file_path = Path(file_name)
        if not file_path.is_absolute():
            results_dir = _CURRENT_DIR / 'Results'
            results_dir.mkdir(exist_ok=True)
            file_path = results_dir / file_path
        if file_path.suffix.lower() != '.xlsx':
            file_path = file_path.with_suffix('.xlsx')

    workbook = Workbook()

    simulation_sheet = workbook.active
    simulation_sheet.title = 'Simulation'
    _build_simulation_sheet(simulation_sheet, pv_project, RESULTS)

    energy_sheet = workbook.create_sheet('Hourly energy')
    _build_energy_sheet(energy_sheet, energy_hourly)

    # if no list is passed, the warnings collected by WarningLog since the start of the run are used
    if warning_messages is None:
        warning_messages = WarningLog.get_messages()

    warnings_sheet = workbook.create_sheet('Warnings')
    _build_warnings_sheet(warnings_sheet, warning_messages)

    try:
        workbook.save(file_path)
    except PermissionError:
        raise PermissionError(f'Could not write "{file_path}". The file is probably already open in Excel, '
                              f'please close it and run the simulation again')

    print(f'\nResults exported to: {file_path}')

    return file_path
