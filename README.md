# PVyield

**Multi-field PV energy yield simulation engine, built on [pvlib](https://github.com/pvlib/pvlib-python).**

PVyield simulates the annual energy yield of utility-scale photovoltaic plants with the goal of reproducing the results of PVsyst as closely as possible, using only open-source tools. It supports projects made of **multiple PV fields**, each with its own modules, inverters, mounting structure, MV transformer and curtailment, all feeding a common HV connection.

> PVyield is an independent project and is not affiliated with or endorsed by PVsyst SA.

## Features

- Multi-field projects (one `PV_Field` per field, aggregated into a `PV_Project`)
- Single-axis trackers (with backtracking) and fixed-tilt structures, 1P/2P/…/1L/2L/… configurations
- Bifacial modules
- Module and inverter data imported directly from **PVsyst `.PAN` and `.OND` files**
- Weather data from **PVGIS** (free, automatic), **Solargis TMY** files, **Solcast API** and **Meteonorm API**
- Horizon shading from PVGIS
- Full loss chain: IAM, soiling, mismatch, DC/AC wiring, inverter, MV/HV transformers, curtailment
- MV cable sizing module (`MV_cable_sizing.py`)
- Excel report with results and the list of warnings raised during the simulation

## Repository layout

| Path | Content |
|---|---|
| `PVyield_launcher.py` | **Start here**: project definition and simulation run |
| `PVyield.py` | Simulation engine |
| `PVyield_classes.py` | Input/output data classes with validation |
| `PVyield_report.py` | Excel report generator |
| `MV_cable_sizing.py` | MV connection cable sizing |
| `classic/` | Previous, dictionary-based version of the engine |
| `PAN/`, `OND/` | Put your module (`.PAN`) and inverter (`.OND`) files here |
| `TMY/` | Put your Solargis TMY files here (only for the `SolarGIS TMY` source) |
| `Results/` | Output reports are written here |

## Installation

```bash
git clone https://github.com/tonival94/PVyield.git
cd PVyield
pip install -r requirements.txt
```

## Usage

1. Copy the `.PAN` file of your module into `PAN/` and the `.OND` file of your inverter into `OND/` (manufacturers usually provide them on request or on their website). Alternatively set `PAN_import = False` / `OND_import = False` and enter the parameters manually.
2. Open `PVyield_launcher.py` and set location, weather source, fields, structures, components and losses. Every parameter is documented inline with its unit.
3. Run it:

```bash
python PVyield_launcher.py
```

The Excel report is saved in `Results/`.

## Inputs

A project is described in `PVyield_launcher.py` by a set of input classes (defined in `PVyield_classes.py`). Each `PV_Field` groups the classes describing one field; `PV_Project` groups the fields together with the parameters shared by the whole plant. Optional inputs left to `None` fall back to the default shown, or are inferred from the other inputs; a warning is raised and listed in the Excel report whenever this happens.

```
PV_Project
├── Coordinates, Tmy
├── PV_Field (one or more)
│   ├── Power, Mms, Ground
│   ├── Module, Inverter
│   ├── MV_Trafo, Field_Curtailment, Loss
│   └── MV_Connection
├── HV_Trafo, HV_Connection (optional)
└── Project_Curtailment (optional)
```

<details>
<summary><b>Site and weather</b> – <code>Coordinates</code>, <code>Tmy</code>, <code>Ground</code></summary>

| Class | Input | Unit | Description |
|---|---|---|---|
| `Coordinates` | `latitude` | deg | Decimal degrees, positive North [-90, 90] |
| | `longitude` | deg | Decimal degrees, positive East [-180, 180] |
| `Tmy` | `TMY_source` | – | `PVGIS`, `SolarGIS TMY`, `Solcast API`, `Meteonorm API` (not case sensitive) |
| | `api_key` | – | API key for subscription-based sources |
| | `TMY_file_name` | – | Solargis TMY file in the `TMY/` folder (only for `SolarGIS TMY`) |
| `Ground` | `albedo` | – | Ground albedo |
| | `slope_tilt` | % | Ground slope |
| | `slope_azimuth` | deg | Downhill direction of the slope (0 = N, 90 = E, 180 = S, 270 = W) |

</details>

<details>
<summary><b>Field size and mounting structure</b> – <code>Power</code>, <code>Mms</code></summary>

| Class | Input | Unit | Description |
|---|---|---|---|
| `Power` | `p_dc_peak` | MWp | Field DC peak power, rounded down to a whole number of strings |
| | `p_ac_nom` | MVA | Field nominal inverter power, rounded up to a whole number of inverters |
| `Mms` | `mount_type` | – | `Single-axis tracker` or `Fixed tilt` |
| | `mod_configuration` | – | `1P`, `2P`, `3P`… (portrait) or `1L`, `2L`… (landscape) |
| | `modules_per_string` | – | Modules per string |
| | `pitch` | m | Row-to-row pitch |
| | `tilt_angle` | deg | Fixed tilt: surface tilt. Tracker: max rotation angle, or a `(min, max)` tuple |
| | `ground_gap` | m | Minimum module-to-ground clearance (for trackers, at max rotation) |
| | `axis_azimuth` | deg | Structure axis azimuth, clockwise from North. Default 180 for trackers, 90 for fixed tilt |
| | `backtrack` | – | Tracker backtracking on/off |

</details>

<details>
<summary><b>PV module</b> – <code>Module</code></summary>

Set `PAN_import = True` and `PAN_file_name` to read every parameter from a PVsyst `.PAN` file in `PAN/`, or `PAN_import = False` to enter them manually:

| Input | Unit | Description |
|---|---|---|
| `manufacturer`, `model` | – | Module identification |
| `Pnom` | Wp | Power at STC |
| `height`, `width` | m | Module dimensions |
| `cells_in_series`, `cells_in_parallel` | – | Cell layout (half-cut: series count divided by 2, parallel = 2) |
| `n_diodes` | – | Bypass diodes |
| `I_L_ref`, `Voc`, `imp_ref`, `vmp_ref` | A, V | Photocurrent, open-circuit voltage and MPP at STC |
| `alpha_sc` | mA/°C | Short-circuit current temperature coefficient (`muIsc` in the PAN) |
| `gamma_ref`, `mu_gamma` | –, 1/K | Diode ideality factor and its temperature coefficient |
| `R_sh_ref`, `R_sh_0`, `R_sh_exp`, `R_s` | Ω, – | Shunt and series resistances of the PVsyst one-diode model (inferred if `None`) |
| `alpha_absorption` | – | Absorptivity |
| `bifaciality` | – | Bifaciality factor |
| `theta_ref`, `iam_ref` | deg, – | Custom IAM profile |
| `V_bypass_diode` | V | Bypass diode voltage drop (default 0.7) |
| `substring_type` | – | `slHalfCells` (half-cut) or `slInLength` (full-cut) |

</details>

<details>
<summary><b>Inverter</b> – <code>Inverter</code></summary>

Set `OND_import = True` and `OND_file_name` to read every parameter from a PVsyst `.OND` file in `OND/` (including efficiency curves at three voltages and temperature derating), or `OND_import = False` to enter them manually:

| Input | Unit | Description |
|---|---|---|
| `manufacturer`, `model` | – | Inverter identification |
| `Pnom_ac`, `P_max_ac` | kW | Nominal and maximum AC power |
| `v_min_mppt`, `v_max_mppt` | V | MPPT voltage window |
| `i_max_mppt` | A | Max current per MPPT |
| `n_string_inputs`, `n_mppt` | – | String inputs and MPPT count |
| `v_out` | V | AC output voltage |
| `inv_euro_eff` | % | European efficiency (manual input only) |
| `T_p_nom`, `T_p_max` | °C | Ambient temperature for nominal (default 40) and maximum (default 30) power |

</details>

<details>
<summary><b>Losses</b> – <code>Loss</code></summary>

| Input | Unit | Description |
|---|---|---|
| `iam_model` | – | `PAN` (IAM profile of the module) or `GlassAR` (Martin-Ruiz, anti-reflective glass) |
| `soiling_loss` | % | Default 2 |
| `rear_shade_factor` | % | Rear side shading from structures (PVsyst standard 5) |
| `thermal_loss_factor` | W/m²K | Module thermal loss factor (default 29) |
| `module_quality` | % | Negative value is a gain (default −0.3) |
| `lid_letid_loss` | % | Default 0.5 |
| `first_year_degrad` | % | Default 0.3 |
| `mismatch_loss_front`, `mismatch_loss_back` | % | Default 0.5 and 0.2 |
| `dc_ohmic_loss_STC` | % | DC wiring loss at STC (default 1.5) |
| `auxiliaries_loss` | % | Inferred from the mounting type if `None` |
| `pv_availability` | % | Default 99 for trackers, 99.5 for fixed tilt |
| `ac_lv_cable_len`, `ac_lv_cable_sect`, `ac_lv_cable_material` | m, mm², – | LV AC cable from inverter to transformer (`Cu` or `Alu`) |
| `ac_mv_cable_len`, `ac_mv_cable_sect`, `ac_mv_cable_material` | m, mm², – | MV AC cable from transformer to collection cabin (`Cu` or `Alu`) |

</details>

<details>
<summary><b>Grid connection and curtailment</b> – <code>MV_Trafo</code>, <code>MV_Connection</code>, <code>HV_Trafo</code>, <code>HV_Connection</code>, <code>Field_Curtailment</code>, <code>Project_Curtailment</code></summary>

| Class | Input | Unit | Description |
|---|---|---|---|
| `MV_Trafo` | `v_out` | kV | LV/MV transformer output voltage (typically 15, 20 or 30) |
| | `p_max` | MVA | Rated power (adapted to the field if `None`) |
| | `iron_loss`, `copper_loss` | % | Default 0.1 and 1.0 |
| `MV_Connection` | `line_len` | km | MV line from the collection cabin to the HV transformer or PoC (`None` = no line) |
| | `material` | – | `Cu` or `Alu` |
| | `destination_field` | – | Name of the field the line connects to, or `grid` / `PoC`: fields can be chained |
| | `auto_sizing` | – | Automatic cable sizing (number of lines and section) |
| | `line_nr`, `line_sect` | –, mm² | Number of lines and section, if not auto-sized |
| `HV_Trafo` | `p_nom` | MVA | MV/HV transformer rated power |
| | `iron_loss`, `copper_loss` | % | Default 0.4 and 0.2 (0.04 above 50 MVA) |
| `HV_Connection` | `line_len` | km | HV line from the HV transformer to the PoC (`None` = no line) |
| `Field_Curtailment` | `grid_curtail_point` | – | `None`, `Inverter out` or `Field MV cabin` |
| | `grid_curtail_power` | MVA | Power limit at that point |
| `Project_Curtailment` | `grid_curtail_point` | – | `None`, `Final MV line in`, `HV trafo in`, `HV trafo out` or `HV line out` |
| | `grid_curtail_power` | MVA | Power limit at that point |
| `PV_Project` | `grid_availability` | % | Default 99.7 |

</details>

## Weather data sources

| `TMY_source` | Requirements |
|---|---|
| `PVGIS` | None: downloaded automatically from the project coordinates |
| `SolarGIS TMY` | A Solargis TMY CSV file in `TMY/` |
| `Solcast API` | Solcast subscription and API key |
| `Meteonorm API` | Meteonorm subscription and API key |

## License

PVyield is released under the **GNU Affero General Public License v3.0** – see [LICENSE](LICENSE).
You can use it freely, including for internal commercial work. If you distribute a modified version or offer it to others as a network service, you must release your source code under the same license.
For commercial licensing under different terms, contact the author.

Copyright © 2024-2026 Antonio Vallone
