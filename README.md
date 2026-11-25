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
