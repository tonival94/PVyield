# -*- coding: utf-8 -*-
"""
MV cable sizing

Python conversion of the Excel VBA module 'Get_Cable_Size_Opt' (MV_Connection_sizing___Losses_v8.xlsm).
Computes the optimal MV cable configuration (number of triads x cross section) of a MV connection line,
by checking both the ampacity constraint (derated cable current carrying capability) and the energy loss
constraint (Monte Carlo simulation of the yearly energy loss along the line).

@author: Antonio Vallone (-T O N I-)

Libraries requirements:
> numpy == 2.4.0

Unit convention:
    - all the arguments exposed to the user keep the engineering units of the original VBA UDF
      ([MW], [MVA], [km], [kV], [mm2], [%], ...), so that the module can be used as a drop-in replacement
    - the cable database keeps the datasheet units ([Ohm/km], [uF/km], [mH/km], [A]), so that it stays
      directly comparable with the manufacturer catalogue, and GetCableProp returns those same units
    - every physical calculation is performed in SI units ([W], [V], [A], [m], [Ohm/m], [H/m]);
      the conversion is done inside each function, right before the calculation
    - two quantities are deliberately NOT converted: the ground temperature stays in [degC] because all the
      derating tables are tabulated in [degC], and the equivalent hours stay in [h] because they only appear
      as a ratio (energy lost / energy produced), where the time unit cancels out

# improvements for next versions:
# 1 - GetEnergyLoss models the energy loss as dV*I with dV = I*(R*cos(phi) + X*sin(phi)). The reactive term
#     X*sin(phi) is a voltage drop, not a power dissipation, so the loss is overestimated. Kept as per VBA
# 2 - add the sizing of DC systems, through a 'System' argument switching between 'ac' and 'dc'. On hold until
#     the DC cable datasheets are available. In DC there is no cos(phi), no reactance and no square root of 3,
#     so the formulas become dV = 2*R*L*I (pole to pole, both conductors in series along the path) and
#     Ploss = 2*N*R*L*I*I, with I = P/V and one circuit made of two conductors instead of a triad of three.
#     Open points to be settled once the datasheets are in:
#       - voltage level: the database below only covers 15 to 36 kV AC XLPE cables, so a separate table is
#         needed for LV DC (1500 V field collection, DC coupled BESS)
#       - ampacity: in DC there is no skin effect, no proximity effect and no dielectric loss, so the real
#         ampacity is higher than the AC one tabulated here. Either keep the AC values (conservative) or apply
#         an uplift factor, which must come from the datasheets and not be guessed
#       - grouping derating: k4 and k5 are tabulated for groups of three cables, while a DC circuit has two.
#         Reusing them as they are is conservative, k1, k2 and k3 stay valid as they only depend on the soil
"""

import math
import warnings
import numpy as np
from typing import List


# ---------------------------------------------------------------------------
#   VALIDATION LISTS
# ---------------------------------------------------------------------------

_VALID_CABLE_CORE     = ['mono', 'multi']       # cable construction: one core per phase, or three cores in one cable
_VALID_CABLE_LAY      = ['trefoil', 'flat']      # cable laydown method inside the trench
_VOLTAGE_CLASSES      = [17.5, 24.0, 33.0, 36.0]     # [kV] rated voltage classes covered by the cable database
_VALID_CABLE_MATERIAL = ['alu', 'cu', 'aluminum', 'copper']  # allowed cable material
_VALID_MV_CABLE_SECT  = [50, 70, 95, 120, 150, 185, 240, 300, 400, 500, 630]  # [mm2] allowed cable section
_EMPTY                  = object()

# ---------------------------------------------------------------------------
#   CABLE DATABASE
# ---------------------------------------------------------------------------
# Based on the ELSEWEDY XLPE unarmoured cable catalogue.
# 25 mm2, 35 mm2 and 800 mm2 are not included: 25/35 are too small to be of any use on a MV connection line
# and 800 is not a standard commercial size in EU. For multi core cables 400/500/630 are not available either.
#
# Single core (mono) database columns:
#   0 - Cross section area                                  [mm2]
#   1 - Max conductor resistance @ 90 degC                   [Ohm/km]
#   2 - Capacitance 'C' @ 50 Hz                              [uF/km]
#   3 - Inductance 'L' @ 50 Hz, trefoil laydown              [mH/km]
#   4 - Inductance 'L' @ 50 Hz, flat laydown                 [mH/km]
#   5 - Ampacity, underground, flat, directly buried         [A] (1 cable per phase, flat, touching 'OOO')
#   6 - Ampacity, underground, trefoil, directly buried      [A] (1 cable per phase, trefoil, touching)
#   7 - Ampacity, underground, trefoil, in duct              [A] (1 cable per phase, trefoil, touching)

_DB_CU_15_MONO = np.array([[50,  0.4937, 0.2104, 0.3932, 0.4394, 238, 245, 183],
                           [70,  0.3422, 0.2421, 0.3696, 0.4158, 292, 295, 227],
                           [95,  0.2467, 0.2620, 0.3568, 0.4030, 347, 351, 272],
                           [120, 0.1959, 0.2851, 0.3456, 0.3918, 395, 398, 314],
                           [150, 0.1592, 0.3081, 0.3348, 0.3810, 443, 445, 355],
                           [185, 0.1278, 0.3360, 0.3223, 0.3685, 501, 502, 412],
                           [240, 0.0981, 0.3752, 0.3104, 0.3567, 580, 578, 478],
                           [300, 0.0792, 0.4177, 0.2994, 0.3457, 658, 648, 547],
                           [400, 0.0633, 0.4600, 0.2908, 0.3370, 746, 726, 623],
                           [500, 0.0511, 0.5153, 0.2814, 0.3276, 840, 817, 711],
                           [630, 0.0417, 0.5713, 0.2739, 0.3201, 940, 906, 806]])

_DB_ALU_15_MONO = np.array([[50,  0.8220, 0.2104, 0.3932, 0.4394, 188, 189, 146],
                            [70,  0.5682, 0.2371, 0.3732, 0.4194, 228, 229, 181],
                            [95,  0.4106, 0.2620, 0.3568, 0.4030, 275, 278, 216],
                            [120, 0.3249, 0.2851, 0.3456, 0.3918, 312, 320, 249],
                            [150, 0.2648, 0.3081, 0.3348, 0.3810, 350, 354, 280],
                            [185, 0.2111, 0.3360, 0.3223, 0.3685, 396, 398, 319],
                            [240, 0.1615, 0.3736, 0.3109, 0.3572, 460, 459, 377],
                            [300, 0.1298, 0.4128, 0.3001, 0.3463, 517, 515, 429],
                            [400, 0.1019, 0.4568, 0.2898, 0.3360, 591, 584, 497],
                            [500, 0.0804, 0.5137, 0.2826, 0.3288, 675, 667, 575],
                            [630, 0.0639, 0.5705, 0.2746, 0.3209, 768, 778, 660]])

_DB_CU_20_MONO = np.array([[50,  0.4937, 0.1812, 0.4121, 0.4583, 238, 245, 183],
                           [70,  0.3421, 0.2072, 0.3854, 0.4317, 292, 297, 227],
                           [95,  0.2467, 0.2235, 0.3734, 0.4196, 347, 353, 275],
                           [120, 0.1959, 0.2424, 0.3599, 0.4061, 397, 402, 318],
                           [150, 0.1591, 0.2612, 0.3496, 0.3959, 443, 448, 358],
                           [185, 0.1277, 0.2839, 0.3351, 0.3813, 501, 505, 416],
                           [240, 0.0980, 0.3159, 0.3222, 0.3684, 580, 579, 481],
                           [300, 0.0790, 0.3505, 0.3114, 0.3577, 662, 650, 550],
                           [400, 0.0631, 0.3850, 0.3020, 0.3482, 750, 729, 625],
                           [500, 0.0508, 0.4300, 0.2916, 0.3378, 844, 819, 711],
                           [630, 0.0414, 0.4755, 0.2833, 0.3295, 945, 906, 807]])

_DB_ALU_20_MONO = np.array([[50,  0.8220, 0.1812, 0.4121, 0.4583, 188, 189, 148],
                            [70,  0.5682, 0.2031, 0.3892, 0.4354, 228, 232, 181],
                            [95,  0.4106, 0.2235, 0.3734, 0.4196, 275, 278, 222],
                            [120, 0.3248, 0.2424, 0.3599, 0.4061, 312, 320, 254],
                            [150, 0.2647, 0.2612, 0.3496, 0.3959, 350, 354, 286],
                            [185, 0.2111, 0.2839, 0.3351, 0.3813, 396, 405, 322],
                            [240, 0.1614, 0.3146, 0.3227, 0.3690, 460, 465, 381],
                            [300, 0.1297, 0.3465, 0.3122, 0.3584, 517, 520, 435],
                            [400, 0.1017, 0.3823, 0.3010, 0.3472, 591, 600, 505],
                            [500, 0.0802, 0.4286, 0.2928, 0.3390, 675, 675, 580],
                            [630, 0.0637, 0.4749, 0.2840, 0.3303, 768, 790, 675]])

_DB_CU_30_MONO = np.array([[50,  0.4937, 0.1412, 0.4508, 0.4970, 238, 245, 183],
                           [70,  0.3421, 0.1596, 0.4215, 0.4677, 292, 297, 230],
                           [95,  0.2466, 0.1711, 0.4078, 0.4540, 347, 356, 279],
                           [120, 0.1958, 0.1844, 0.3927, 0.4390, 397, 404, 322],
                           [150, 0.1590, 0.1976, 0.3809, 0.4272, 447, 451, 362],
                           [185, 0.1275, 0.2135, 0.3648, 0.4111, 505, 509, 420],
                           [240, 0.0977, 0.2358, 0.3499, 0.3961, 585, 583, 485],
                           [300, 0.0787, 0.2599, 0.3372, 0.3834, 666, 654, 554],
                           [400, 0.0627, 0.2839, 0.3261, 0.3723, 754, 731, 629],
                           [500, 0.0503, 0.3151, 0.3138, 0.3601, 848, 821, 713],
                           [630, 0.0408, 0.3467, 0.3039, 0.3501, 949, 908, 810]])

_DB_ALU_30_MONO = np.array([[50,  0.8220, 0.1412, 0.4508, 0.4970, 188, 189, 152],
                            [70,  0.5682, 0.1567, 0.4256, 0.4719, 228, 232, 189],
                            [95,  0.4106, 0.1711, 0.4078, 0.4540, 279, 278, 222],
                            [120, 0.3248, 0.1844, 0.3927, 0.4390, 316, 320, 259],
                            [150, 0.2646, 0.1976, 0.3809, 0.4272, 354, 354, 290],
                            [185, 0.2110, 0.2135, 0.3648, 0.4111, 400, 405, 322],
                            [240, 0.1612, 0.2349, 0.3505, 0.3967, 460, 468, 386],
                            [300, 0.1295, 0.2571, 0.3382, 0.3844, 517, 526, 440],
                            [400, 0.1015, 0.2820, 0.3252, 0.3714, 591, 605, 510],
                            [500, 0.0799, 0.3142, 0.3150, 0.3612, 679, 684, 587],
                            [630, 0.0632, 0.3462, 0.3046, 0.3508, 772, 794, 680]])

# Multi core database columns:
#   0 - Cross section area                                  [mm2]
#   1 - Max conductor resistance @ 90 degC                   [Ohm/km]
#   2 - Capacitance 'C' @ 50 Hz                              [uF/km]
#   3 - Inductance 'L' @ 50 Hz                               [mH/km]
#   4 - Ampacity, underground, directly buried               [A]
#   5 - Ampacity, underground, in duct                       [A]

_DB_CU_15_MULTI = np.array([[50,  0.4938, 0.2104, 0.3715, 234, 181],
                            [70,  0.3422, 0.2421, 0.3482, 286, 223],
                            [95,  0.2468, 0.2620, 0.3366, 343, 268],
                            [120, 0.1960, 0.2851, 0.3252, 390, 307],
                            [150, 0.1594, 0.3081, 0.3156, 438, 346],
                            [185, 0.1280, 0.3360, 0.3030, 496, 394],
                            [240, 0.0984, 0.3752, 0.2917, 575, 461],
                            [300, 0.0796, 0.4177, 0.2824, 657, 530]])

_DB_ALU_15_MULTI = np.array([[50,  0.8220, 0.2104, 0.3715, 181, 140],
                             [70,  0.5683, 0.2371, 0.3514, 222, 173],
                             [95,  0.4107, 0.2620, 0.3366, 266, 208],
                             [120, 0.3249, 0.2851, 0.3252, 303, 238],
                             [150, 0.2649, 0.3081, 0.3156, 340, 269],
                             [185, 0.2113, 0.3360, 0.3030, 386, 307],
                             [240, 0.1617, 0.3736, 0.2921, 449, 360],
                             [300, 0.1300, 0.4128, 0.2829, 509, 411]])

_DB_CU_20_MULTI = np.array([[50,  0.4937, 0.1812, 0.3917, 234, 184],
                            [70,  0.3422, 0.2072, 0.3667, 287, 226],
                            [95,  0.2467, 0.2235, 0.3541, 343, 272],
                            [120, 0.1960, 0.2424, 0.3417, 391, 311],
                            [150, 0.1592, 0.2612, 0.3312, 439, 351],
                            [185, 0.1279, 0.2839, 0.3177, 496, 399],
                            [240, 0.0982, 0.3159, 0.3052, 576, 466],
                            [300, 0.0794, 0.3505, 0.2948, 657, 536]])

_DB_ALU_20_MULTI = np.array([[50,  0.8220, 0.1812, 0.3917, 181, 142],
                             [70,  0.5682, 0.2031, 0.3702, 223, 176],
                             [95,  0.4107, 0.2235, 0.3541, 266, 211],
                             [120, 0.3249, 0.2424, 0.3417, 304, 242],
                             [150, 0.2648, 0.2612, 0.3312, 341, 272],
                             [185, 0.2112, 0.2839, 0.3177, 386, 311],
                             [240, 0.1615, 0.3146, 0.3057, 450, 364],
                             [300, 0.1299, 0.3465, 0.2954, 509, 415]])

_DB_CU_30_MULTI = np.array([[50,  0.4937, 0.1412, 0.4333, 236, 191],
                            [70,  0.3421, 0.1596, 0.4052, 287, 233],
                            [95,  0.2466, 0.1711, 0.3909, 344, 280],
                            [120, 0.1958, 0.1844, 0.3766, 391, 319],
                            [150, 0.1591, 0.1976, 0.3645, 439, 360],
                            [185, 0.1276, 0.2135, 0.3492, 497, 409],
                            [240, 0.0979, 0.2358, 0.3344, 577, 477],
                            [300, 0.0790, 0.2599, 0.3219, 657, 547]])

_DB_ALU_30_MULTI = np.array([[50,  0.8220, 0.1412, 0.4333, 183, 148],
                             [70,  0.5682, 0.1567, 0.4091, 223, 181],
                             [95,  0.4106, 0.1711, 0.3909, 266, 217],
                             [120, 0.3248, 0.1844, 0.3766, 304, 248],
                             [150, 0.2647, 0.1976, 0.3645, 341, 279],
                             [185, 0.2110, 0.2135, 0.3492, 387, 318],
                             [240, 0.1613, 0.2349, 0.3350, 450, 372],
                             [300, 0.1297, 0.2571, 0.3228, 509, 424]])

# Cable database lookup table, keyed by (material, voltage class [kV], core).
# The 36 kV entries are a copy of the 30 kV ones: given the lack of published data, 36 kV copper cables are
# assumed to have the same electrical properties as the 30 kV ones. Aluminium is NOT available at 36 kV.
_CABLE_DATABASE = {('cu',  17.5, 'mono'):  _DB_CU_15_MONO,
                   ('alu', 17.5, 'mono'):  _DB_ALU_15_MONO,
                   ('cu',  24.0, 'mono'):  _DB_CU_20_MONO,
                   ('alu', 24.0, 'mono'):  _DB_ALU_20_MONO,
                   ('cu',  33.0, 'mono'):  _DB_CU_30_MONO,
                   ('alu', 33.0, 'mono'):  _DB_ALU_30_MONO,
                   ('cu',  36.0, 'mono'):  _DB_CU_30_MONO,
                   ('cu',  17.5, 'multi'): _DB_CU_15_MULTI,
                   ('alu', 17.5, 'multi'): _DB_ALU_15_MULTI,
                   ('cu',  24.0, 'multi'): _DB_CU_20_MULTI,
                   ('alu', 24.0, 'multi'): _DB_ALU_20_MULTI,
                   ('cu',  33.0, 'multi'): _DB_CU_30_MULTI,
                   ('alu', 33.0, 'multi'): _DB_ALU_30_MULTI,
                   ('cu',  36.0, 'multi'): _DB_CU_30_MULTI,
                   }


# ---------------------------------------------------------------------------
#   DERATING FACTOR TABLES
# ---------------------------------------------------------------------------
# All the derating tables of the original VBA code were coded as chained Select Case blocks performing a
# linear interpolation between two consecutive tabulated points. They are stored here as breakpoint arrays
# and interpolated with np.interp, which gives numerically identical results in a fraction of the code.

# k1 - ground temperature derating factor
_K1_TEMPERATURE = np.array([15, 20, 25, 30, 35, 40, 45, 50, 55])                    # [degC] ground temperature
_K1_FACTOR = np.array([1.04, 1.00, 0.96, 0.93, 0.89, 0.85, 0.80, 0.76, 0.71])       # [-] derating factor

# k2 - trench depth derating factor. One row per installation case
_K2_DEPTH = np.array([0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0])          # [m] trench depth
_K2_FACTOR = {('buried', 'mono', 'small'): np.array([1.00, 0.98, 0.96, 0.94, 0.92, 0.91, 0.90, 0.89, 0.88, 0.87]),
              ('buried', 'mono', 'large'): np.array([1.00, 0.98, 0.94, 0.92, 0.90, 0.88, 0.86, 0.85, 0.83, 0.81]),
              ('buried', 'multi', 'any'):  np.array([1.00, 0.99, 0.96, 0.94, 0.92, 0.91, 0.90, 0.89, 0.88, 0.87]),
              ('duct',   'mono', 'small'): np.array([1.00, 0.98, 0.96, 0.94, 0.92, 0.91, 0.90, 0.89, 0.88, 0.87]),
              ('duct',   'mono', 'large'): np.array([1.00, 0.98, 0.95, 0.92, 0.90, 0.89, 0.88, 0.87, 0.85, 0.84]),
              ('duct',   'multi', 'any'):  np.array([1.00, 0.99, 0.97, 0.96, 0.94, 0.93, 0.92, 0.91, 0.90, 0.89]),
              }
_K2_AREA_THRESHOLD = 185        # [mm2] above this cross section the 'large' derating table applies

# k3 - soil thermal resistivity derating factor
_K3_RESISTIVITY = np.array([0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0])                # [K*m/W] soil resistivity
_K3_FACTOR = np.array([1.10, 1.05, 1.00, 0.92, 0.83, 0.73, 0.66, 0.60])             # [-] derating factor

# k4 - grouping derating factor for directly buried cables. One matrix per (laydown, core) combination,
# with one row per number of triads (index 0 -> 1 triad ... index 5 -> 6 or more triads)
_K4_SPACING = np.array([0.0, 0.15, 0.30])                                           # [m] spacing between triads
_K4_FACTOR = {('trefoil', 'mono'): np.array([[1.00, 1.00, 1.00],
                                             [0.77, 0.82, 0.88],
                                             [0.66, 0.73, 0.80],
                                             [0.60, 0.68, 0.74],
                                             [0.56, 0.64, 0.72],
                                             [0.53, 0.61, 0.70]]),
              ('flat', 'mono'):    np.array([[1.00, 1.00, 1.00],
                                             [0.80, 0.85, 0.91],
                                             [0.69, 0.76, 0.83],
                                             [0.63, 0.71, 0.77],
                                             [0.59, 0.67, 0.75],
                                             [0.57, 0.64, 0.73]]),
              ('trefoil', 'multi'): np.array([[1.00, 1.00, 1.00],
                                              [0.81, 0.87, 0.91],
                                              [0.69, 0.76, 0.82],
                                              [0.62, 0.72, 0.77],
                                              [0.58, 0.66, 0.73],
                                              [0.54, 0.63, 0.70]]),
              ('flat', 'multi'):    np.array([[1.00, 1.00, 1.00],
                                              [0.81, 0.87, 0.91],
                                              [0.70, 0.78, 0.84],
                                              [0.63, 0.74, 0.81],
                                              [0.60, 0.70, 0.78],
                                              [0.56, 0.67, 0.76]]),
              }
_K4_MAX_TRIADS = 6              # [-] the grouping tables are tabulated up to 6 triads laid at the same depth

# k5 - grouping derating factor for cables laid in a duct
_K5_SIX_TRIADS = 0.51           # [-] derating factor for exactly 6 triads laid in a duct

_MAX_TRENCH_DEPTH = 3.0         # [m] deepest trench covered by the k2 derating table




# ---------------------------------------------------------------------------
#   MAIN FUNCTION
# ---------------------------------------------------------------------------

def GetCableSizeOpt(Pdc,
                    Pac,
                    EqHours,
                    Lenght,
                    V,
                    Material,
                    MaxEloss=1.5,
                    Core='mono',
                    Lay='trefoil',
                    Cosphi=0.9,
                    SoilR=1.5,
                    GroundT=30,
                    Depth=0.8,
                    Spacing=0.15,
                    Duct=False,
                    MaxLoading=0.8,
                    BetaA=0.53,
                    BetaB=0.57,
                    f50=True,
                    f70=True,
                    f95=True,
                    f120=True,
                    f150=True,
                    f185=True,
                    f240=True,
                    f300=True,
                    f400=True,
                    f500=True,
                    f630=True,
                    iterations=3,
                    seed=None,
                    freq=50,
                    ):
    """
    Computes the optimal MV cable configuration of a connection line, i.e. the number of cable triads and the
    cross section that carry the plant power to the point of connection at the lowest cost.

    For every allowed cross section, two independent constraints are evaluated:
        - ampacity: how many triads are needed so that the derated ampacity of the cable covers the max line
          current. Up to 6 triads can be laid at the same depth, further triads are pushed one metre deeper
        - energy loss: how many triads are needed so that the yearly energy loss along the line stays below
          the admissible one. Since the loss scales with the inverse of the number of triads, it is computed
          once on a trial configuration and then rescaled
    The number of triads of each cross section is the max between the two, and the optimal configuration is
    the smallest cross section that still requires the minimum number of triads.

    Args:
        Pdc (float):          [MWp] DC power. Used only for PV plants. For BESS, DC power equals AC power
        Pac (float):          [MVA] AC power at the injection point
        EqHours (float):      [kWh/kWp or h] plant yield (PV) or equivalent full load hours (BESS)
        Lenght (float):       [km] connection line length
        V (float):            [kV] connection line rated phase-to-phase voltage
        Material (str):       [-] cable material, either 'cu'/'copper' or 'alu'/'aluminum'. Not case sensitive
        MaxEloss (float):     [%] max admissible energy loss along the connection line. Default 1.5
        Core (str):           [-] cable construction, either 'mono' or 'multi'. Default 'mono'
        Lay (str):            [-] cable laydown method, either 'trefoil' or 'flat'. Default 'trefoil'
        Cosphi (float):       [-] power factor along the connection line. Default 0.9
        SoilR (float):        [K*m/W] soil thermal resistivity. Default 1.5
        GroundT (float):      [degC] ground temperature. Default 30
        Depth (float):        [m] trench depth. Default 0.8
        Spacing (float):      [m] spacing between two adjacent cable triads. Default 0.15
        Duct (bool):          [-] True if the cables are laid inside a duct. Default True
        MaxLoading (float):   [-] max loading allowed on the cable, as a safety margin against the short
                                  circuit current, which the derating factors do not account for. Default 0.8
        BetaA (float):        [-] alpha parameter of the Beta distribution of the plant power output. Default 0.53
        BetaB (float):        [-] beta parameter of the Beta distribution of the plant power output. Default 0.57
        f50 (bool):           [-] flag: allow the use of a 50 mm2 cable. Default True
        f70 (bool):           [-] flag: allow the use of a 70 mm2 cable. Default True
        f95 (bool):           [-] flag: allow the use of a 95 mm2 cable. Default True
        f120 (bool):          [-] flag: allow the use of a 120 mm2 cable. Default True
        f150 (bool):          [-] flag: allow the use of a 150 mm2 cable. Default True
        f185 (bool):          [-] flag: allow the use of a 185 mm2 cable. Default True
        f240 (bool):          [-] flag: allow the use of a 240 mm2 cable. Default True
        f300 (bool):          [-] flag: allow the use of a 300 mm2 cable. Default True
        f400 (bool):          [-] flag: allow the use of a 400 mm2 cable, single core only. Default True
        f500 (bool):          [-] flag: allow the use of a 500 mm2 cable, single core only. Default True
        f630 (bool):          [-] flag: allow the use of a 630 mm2 cable, single core only. Default True
        iterations (int):     [-] number of Monte Carlo runs used to smooth the energy loss estimate. Default 3
        seed (int or None):   [-] seed of the random number generator. Default None (non reproducible run)
        freq (float):         [Hz] grid frequency. Default 50

    Returns:
        n_triads (int):       [-] optimal number of cable triads of the connection line
        area (float):         [mm2] optimal cable cross section of the connection line
    """

    # 1 - INPUT CHECK
    (Pdc, Pac, EqHours, Lenght, V, Material, MaxEloss, Core, Lay, Cosphi, SoilR, GroundT,
     Depth, Spacing, Duct, MaxLoading, BetaA, BetaB, flags, iterations, seed, freq) = \
        GetCableSizeOpt_input_check(Pdc, Pac, EqHours, Lenght, V, Material, MaxEloss, Core, Lay, Cosphi,
                                    SoilR, GroundT, Depth, Spacing, Duct, MaxLoading, BetaA, BetaB,
                                    f50, f70, f95, f120, f150, f185, f240, f300, f400, f500, f630,
                                    iterations, seed, freq)

    # 2 - UNIT CONVERSION
    # the max admissible loss is the only argument converted here: it is the only one used directly by this
    # function. Every other quantity is converted to SI inside the subfunction that actually uses it
    max_eloss = MaxEloss / 100                                  # [-] max admissible energy loss, as a fraction in the range [0, 1]
    p_ac_si = Pac * 1E6                                         # [VA] AC power at the injection point
    v_si = V * 1E3                                              # [V] phase-to-phase voltage of the connection line

    # 3 - LOAD THE CABLE DATABASE
    database = GetCableDatabase(Material, V, Core, *flags)      # [various] cable properties of every allowed cross section
    areas = database[:, 0]                                      # [mm2] cross sections to be evaluated, sorted from the smallest to the largest
    n_areas = len(areas)                                        # [-] number of cross sections to be evaluated

    # 4 - COMPUTE NtA, the number of cable triads required by the cable ampacity
    nt_ampacity = np.zeros(n_areas)                             # [-] preallocate the number of triads allowed by the ampacity

    for i, area in enumerate(areas):

        i_line = p_ac_si / (math.sqrt(3) * v_si)                                    # [A] max line current flowing from the plant to the point of connection
        ampacity = GetCableProp('A', Material, V, area, Core, Lay, Duct)            # [A] nominal ampacity of one cable, before any derating
        current_depth = Depth                                                       # [m] depth of the layer being filled with triads
        n_triads = 0                                                                # [-] triads accumulated over all the layers filled so far
        j = 1                                                                       # [-] number of triads assumed to be laid at the current depth

        while True:

            # derate the nominal ampacity for the current installation conditions, then apply the safety margin
            derating = GetDerating(j, area, Core, Lay, SoilR, GroundT, current_depth, Spacing, Duct)    # [-] cumulative derating factor
            corrected_ampacity = derating * MaxLoading * ampacity                                       # [A] ampacity actually usable by one cable

            # number of triads needed to carry the outstanding line current, rounded up to the next integer
            nt_required = math.ceil(i_line / corrected_ampacity) if corrected_ampacity > 0 else np.inf   # [-] triads required at the current depth

            # the solution is consistent as soon as the triads required are covered by the triads assumed to
            # compute the derating. Since the derating factor decreases monotonically with the number of triads,
            # nt_required can only grow with j, so the first j satisfying this test always satisfies it as an
            # equality. The test is written as '<=' rather than '==' so that the loop cannot be trapped in a
            # deeper layer should a derating table ever be edited into a non monotonic one
            if nt_required <= j:
                n_triads += j                                                       # [-] store the triads of this layer and stop
                break

            # the grouping table stops at 6 triads per layer: if no solution is found, this layer is filled with
            # 6 triads and the outstanding current is pushed to a new layer, one metre deeper
            if j == _K4_MAX_TRIADS:
                i_line = i_line - _K4_MAX_TRIADS * corrected_ampacity               # [A] line current still to be carried by the deeper layers
                n_triads += _K4_MAX_TRIADS                                          # [-] this layer is full
                current_depth = current_depth + 1                                   # [m] move one metre deeper
                j = 0                                                               # [-] restart the triad count on the new layer

                # below 3 m no derating factor is tabulated, so this cross section is declared unfeasible
                if current_depth > _MAX_TRENCH_DEPTH:
                    n_triads = np.inf                                               # [-] flag the cross section as unfeasible
                    break

            j += 1                                                                  # [-] try with one more triad at the current depth

        nt_ampacity[i] = n_triads                                                   # [-] store the result for this cross section

    # 5 - COMPUTE NtL, the number of cable triads required by the max admissible energy loss
    # 6 - COMPUTE Nt, the number of cable triads, and
    # 7 - FIND THE OPTIMAL CABLE.
    # Steps 5 to 7 share the same loop, so that as soon as the best cable is found the remaining, smaller,
    # cross sections are not evaluated at all
    nt_guess = 10                                               # [-] trial number of triads used to compute the energy loss
    optimal_nt = np.inf                                         # [-] best number of triads found so far, intentionally very high
    optimal_area = None                                         # [mm2] cross section associated with the best number of triads

    # scan the cross sections from the largest to the smallest: the number of triads can only grow
    for i in range(n_areas - 1, -1, -1):

        area = areas[i]                                                                                                     # [mm2] cross section being evaluated

        # energy loss of the trial configuration, then rescaled: the loss is inversely proportional to the
        # number of triads, so the triads needed to meet the target are a simple proportion
        loss_guess = GetEnergyLoss(nt_guess, area, Pdc, Pac, EqHours, Lenght, V, Material, Core, Lay,
                                   Cosphi, Duct, BetaA, BetaB, iterations, seed, freq)                                      # [-] energy loss of the trial configuration, as a fraction
        nt_loss = math.ceil(nt_guess * loss_guess / max_eloss)                                                              # [-] triads required to meet the max admissible loss

        # the binding constraint is the most demanding of the two
        nt_total = max(nt_ampacity[i], nt_loss)                                                                             # [-] triads required by this cross section

        # an unfeasible cross section stops the scan: every smaller one would be unfeasible as well
        if not np.isfinite(nt_total):
            break

        # keep shrinking the cross section as long as it does not cost extra triads, otherwise stop
        if nt_total <= optimal_nt:
            optimal_nt = nt_total                                                                                           # [-] new best number of triads
            optimal_area = area                                                                                             # [mm2] new best cross section
        else:
            break

    # no cross section among the allowed ones can carry the power within the given constraints
    if optimal_area is None:
        raise ValueError('MV connection: no feasible cable configuration found. Please allow larger cable '
                         'sections, raise the max admissible energy loss or split the plant on more lines')

    # 8 - OUTPUT THE BEST CONFIGURATION
    return int(optimal_nt), float(optimal_area)


# ---------------------------------------------------------------------------
#   ELECTRICAL CALCULATIONS
# ---------------------------------------------------------------------------

def GetEnergyLoss(Nt,
                  Area,
                  Pdc,
                  Pac,
                  EqHours,
                  Lenght,
                  V,
                  Material,
                  Core,
                  Lay,
                  Cosphi,
                  Duct,
                  BetaA=0.53,
                  BetaB=0.57,
                  iterations=3,
                  seed=None,
                  freq=50,
                  ):
    """
    Computes the energy lost along the connection line over the whole plant production, as a fraction of the
    energy produced, with a Monte Carlo simulation.

    The plant power output over the year is modelled as a Beta distribution: at every hour a power level is
    drawn from it, the resulting phase current and voltage drop are computed, and the corresponding energy
    loss is accumulated until the energy produced matches the plant yield. The whole process is repeated
    'iterations' times and the results are averaged, to smooth the Monte Carlo noise.

    Args:
        Nt (int):             [-] number of cable triads of the connection line
        Area (float):         [mm2] cable cross section
        Pdc (float):          [MWp] DC power. Used only for PV plants. For BESS, DC power equals AC power
        Pac (float):          [MVA] AC power at the injection point
        EqHours (float):      [kWh/kWp or h] plant yield (PV) or equivalent full load hours (BESS)
        Lenght (float):       [km] connection line length
        V (float):            [kV] connection line rated phase-to-phase voltage
        Material (str):       [-] cable material, either 'cu'/'copper' or 'alu'/'aluminum'. Not case sensitive
        Core (str):           [-] cable construction, either 'mono' or 'multi'. Not case sensitive
        Lay (str):            [-] cable laydown method, either 'trefoil' or 'flat'. Not case sensitive
        Cosphi (float):       [-] power factor along the connection line
        Duct (bool):          [-] True if the cables are laid inside a duct, False if directly buried
        BetaA (float):        [-] alpha parameter of the Beta distribution of the plant power output. Default 0.53
        BetaB (float):        [-] beta parameter of the Beta distribution of the plant power output. Default 0.57
        iterations (int):     [-] number of Monte Carlo runs used to smooth the estimate. Default 3
        seed (int or None):   [-] seed of the random number generator. Default None (non reproducible run)
        freq (float):         [Hz] grid frequency. Default 50

    Returns:
        energy_loss (float):  [-] energy lost along the connection line, as a fraction in the range [0, 1]
    """

    # 1 - COMPUTE THE ELECTRICAL PROPERTIES OF THE CABLE
    resistance = GetCableProp('R', Material, V, Area, Core, Lay, Duct)      # [Ohm/km] max conductor resistance at 90 degC
    inductance = GetCableProp('L', Material, V, Area, Core, Lay, Duct)      # [mH/km] cable inductance, according to the laydown method

    # 2 - UNIT CONVERSION TO SI
    r_si = resistance / 1E3                                     # [Ohm/m] conductor resistance
    l_si = inductance * 1E-6                                    # [H/m] cable inductance, from mH/km
    x_si = 2 * np.pi * freq * l_si                              # [Ohm/m] cable inductive reactance
    length_si = Lenght * 1E3                                    # [m] connection line length
    p_ac_si = Pac * 1E6                                         # [VA] AC power at the injection point
    p_dc_si = Pdc * 1E6                                         # [W] DC power of the plant
    v_si = V * 1E3                                              # [V] phase-to-phase voltage of the connection line
    sin_phi = math.sqrt(1 - Cosphi ** 2)                        # [-] reactive factor, paired with the given power factor

    # 3 - PRECOMPUTE THE CONSTANT PART OF THE LOSS
    # the phase current of one triad at full plant power. During the year the actual current is this value
    # scaled by the power level drawn from the Beta distribution
    i_phase_max = p_ac_si / (math.sqrt(3) * v_si * Nt)           # [A] phase current of one triad at full power

    # power lost by the whole line at full plant power. The voltage drop is taken as the real part of the
    # voltage drop phasor, since the imaginary part can be recovered by rephasing
    p_loss_max = 3 * Nt * i_phase_max * (r_si * Cosphi + x_si * sin_phi) * length_si * i_phase_max   # [W] line loss at full power

    # equivalent hours accumulated by one hour at full plant power. The AC power is measured at the injection
    # point while the yield is referred to the DC peak power, hence the ratio
    heq_max = Pac / Pdc                                         # [h] equivalent hours produced by one hour at full power

    # 4 - RUN THE MONTE CARLO SIMULATION
    rng = np.random.default_rng(seed)                           # random number generator, seeded for reproducibility
    loss_mc = 0.0                                               # [-] accumulator of the loss over the Monte Carlo runs

    # first guess of how many hours are needed to reach the plant yield, with a 25% margin plus a fixed offset
    # to make the guess robust also on very short yields
    n_hours_guess = int(np.ceil(EqHours / (heq_max * BetaA / (BetaA + BetaB)) * 1.25)) + 100     # [h] size of the sample block

    for _ in range(iterations):

        # draw one power level per hour. Drawing directly from the Beta distribution is the same as applying
        # the inverse Beta cumulative distribution to a uniform random variable, as the original VBA did
        power_level = rng.beta(BetaA, BetaB, n_hours_guess)                     # [-] plant power output of each hour, as a fraction of the AC power

        # cumulate the equivalent hours produced, hour by hour
        heq_cumulated = np.cumsum(power_level * heq_max)                        # [h] equivalent hours accumulated up to each hour

        # keep drawing until the sample block is long enough to reach the plant yield
        while heq_cumulated[-1] < EqHours:
            power_level = np.concatenate([power_level, rng.beta(BetaA, BetaB, n_hours_guess)])       # [-] extend the sample block
            heq_cumulated = np.cumsum(power_level * heq_max)                                          # [h] recompute the cumulated hours

        # index of the first hour that brings the cumulated production up to the plant yield. That hour is
        # included, exactly as the original While loop did, since the check is performed before adding it
        last_hour = int(np.searchsorted(heq_cumulated, EqHours, side='left'))   # [-] index of the last hour of the simulated year

        # the loss of each hour scales with the square of the power level, so the yearly loss is the loss at
        # full power times the sum of the squared power levels
        energy_lost = p_loss_max * np.sum(power_level[:last_hour + 1] ** 2)     # [Wh] energy lost along the line over the simulated year

        # accumulate the loss of this run, normalised on the energy produced by the plant
        loss_mc += energy_lost / (EqHours * p_dc_si)                            # [-] loss of this run, as a fraction

    # average over the runs to smooth the Monte Carlo noise
    return loss_mc / iterations


def GetVLoss(Nt,
             Area,
             Pac,
             Lenght,
             V,
             Material,
             Core,
             Lay,
             Cosphi,
             Duct,
             freq=50,
             ):
    """
    Computes the phase-to-phase voltage drop along the connection line, at the max plant power.

    Args:
        Nt (int):             [-] number of cable triads of the connection line
        Area (float):         [mm2] cable cross section
        Pac (float):          [MVA] AC power at the injection point
        Lenght (float):       [km] connection line length
        V (float):            [kV] connection line rated phase-to-phase voltage
        Material (str):       [-] cable material, either 'cu'/'copper' or 'alu'/'aluminum'. Not case sensitive
        Core (str):           [-] cable construction, either 'mono' or 'multi'. Not case sensitive
        Lay (str):            [-] cable laydown method, either 'trefoil' or 'flat'. Not case sensitive
        Cosphi (float):       [-] power factor along the connection line
        Duct (bool):          [-] True if the cables are laid inside a duct, False if directly buried
        freq (float):         [Hz] grid frequency. Default 50

    Returns:
        v_loss (float):       [V] phase-to-phase voltage drop along the line, at the max phase current
    """

    # 1 - COMPUTE THE ELECTRICAL PROPERTIES OF THE CABLE
    resistance = GetCableProp('R', Material, V, Area, Core, Lay, Duct)      # [Ohm/km] max conductor resistance at 90 degC
    inductance = GetCableProp('L', Material, V, Area, Core, Lay, Duct)      # [mH/km] cable inductance, according to the laydown method

    # 2 - UNIT CONVERSION TO SI
    r_si = resistance / 1E3                                     # [Ohm/m] conductor resistance
    l_si = inductance * 1E-6                                    # [H/m] cable inductance, from mH/km
    x_si = 2 * np.pi * freq * l_si                              # [Ohm/m] cable inductive reactance
    length_si = Lenght * 1E3                                    # [m] connection line length
    p_ac_si = Pac * 1E6                                         # [VA] AC power at the injection point
    v_si = V * 1E3                                              # [V] phase-to-phase voltage of the connection line
    sin_phi = math.sqrt(1 - Cosphi ** 2)                        # [-] reactive factor, paired with the given power factor

    # 3 - COMPUTE THE VOLTAGE DROP
    i_phase_max = p_ac_si / (math.sqrt(3) * v_si * Nt)          # [A] phase current of one triad at full plant power

    # the voltage drop is taken as the real part of the voltage drop phasor, since the imaginary part can be
    # recovered by rephasing. The square root of 3 turns the phase-to-ground drop into a phase-to-phase one
    return math.sqrt(3) * i_phase_max * (r_si * Cosphi + x_si * sin_phi) * length_si


# ---------------------------------------------------------------------------
#   CABLE DATABASE ACCESS
# ---------------------------------------------------------------------------

def GetCableDatabase(Material,
                     V,
                     Core,
                     f50=True,
                     f70=True,
                     f95=True,
                     f120=True,
                     f150=True,
                     f185=True,
                     f240=True,
                     f300=True,
                     f400=True,
                     f500=True,
                     f630=True,
                     ):
    """
    Returns the cable property database of the requested material, voltage class and construction, restricted
    to the cross sections allowed by the flags.

    Args:
        Material (str):       [-] cable material, either 'cu'/'copper' or 'alu'/'aluminum'. Not case sensitive
        V (float):            [kV] connection line rated phase-to-phase voltage. The database is organised in
                                   voltage classes (17.5, 24, 33, 36 kV) and the first class covering V is used
        Core (str):           [-] cable construction, either 'mono' or 'multi'. Not case sensitive
        f50 (bool):           [-] flag: allow the use of a 50 mm2 cable. Default True
        f70 (bool):           [-] flag: allow the use of a 70 mm2 cable. Default True
        f95 (bool):           [-] flag: allow the use of a 95 mm2 cable. Default True
        f120 (bool):          [-] flag: allow the use of a 120 mm2 cable. Default True
        f150 (bool):          [-] flag: allow the use of a 150 mm2 cable. Default True
        f185 (bool):          [-] flag: allow the use of a 185 mm2 cable. Default True
        f240 (bool):          [-] flag: allow the use of a 240 mm2 cable. Default True
        f300 (bool):          [-] flag: allow the use of a 300 mm2 cable. Default True
        f400 (bool):          [-] flag: allow the use of a 400 mm2 cable, single core only. Default True
        f500 (bool):          [-] flag: allow the use of a 500 mm2 cable, single core only. Default True
        f630 (bool):          [-] flag: allow the use of a 630 mm2 cable, single core only. Default True

    Returns:
        database (np.ndarray): [various] cable properties, one row per allowed cross section, sorted by
                                         growing cross section. Columns for single core cables are:
                                         0 - cross section [mm2]
                                         1 - max conductor resistance @ 90 degC [Ohm/km]
                                         2 - capacitance 'C' @ 50 Hz [uF/km]
                                         3 - inductance 'L' @ 50 Hz, trefoil [mH/km]
                                         4 - inductance 'L' @ 50 Hz, flat [mH/km]
                                         5 - ampacity, flat, directly buried [A]
                                         6 - ampacity, trefoil, directly buried [A]
                                         7 - ampacity, trefoil, in duct [A]
                                         Columns for multi core cables are:
                                         0 - cross section [mm2]
                                         1 - max conductor resistance @ 90 degC [Ohm/km]
                                         2 - capacitance 'C' @ 50 Hz [uF/km]
                                         3 - inductance 'L' @ 50 Hz [mH/km]
                                         4 - ampacity, directly buried [A]
                                         5 - ampacity, in duct [A]
    """

    # 1 - NORMALISE THE STRING ARGUMENTS
    material = _assign_str_from_list(Material, 'MV connection: cable material', _VALID_CABLE_MATERIAL)   # validate the material
    material = {'copper': 'cu', 'aluminum': 'alu'}.get(material, material)                              # normalise the material aliases
    core = _assign_str_from_list(Core, 'MV connection: cable core', _VALID_CABLE_CORE)                  # validate the cable construction

    # 2 - PICK THE VOLTAGE CLASS
    # the database is tabulated by voltage class, so the first class that covers the requested voltage is used
    voltage_class = next((v_class for v_class in _VOLTAGE_CLASSES if V <= v_class), None)                # [kV] voltage class covering V

    # a voltage above the highest class has no cable in the database
    if voltage_class is None:
        raise ValueError(f'MV connection: no cable available at {V} kV. The database covers up to '
                         f'{_VOLTAGE_CLASSES[-1]} kV')

    # 3 - LOAD THE FULL DATABASE
    full_database = _CABLE_DATABASE.get((material, voltage_class, core))                                # [various] every cross section of this cable family

    # aluminium is not manufactured at 36 kV, so this combination has no entry in the database
    if full_database is None:
        raise ValueError(f'MV connection: no {material} {core} core cable available at {V} kV. '
                         f'Please select a different material or lower the connection voltage')

    # 4 - REBUILD THE DATABASE WITH THE ALLOWED CROSS SECTIONS ONLY
    flags = [f50, f70, f95, f120, f150, f185, f240, f300, f400, f500, f630]                             # [-] one flag per commercial cross section
    allowed = [sect for sect, flag in zip(_VALID_MV_CABLE_SECT, flags) if flag]                         # [mm2] cross sections enabled by the user

    # cross sections above 300 mm2 are not manufactured as multi core cables, so they are silently dropped
    unavailable = [sect for sect in allowed if sect not in full_database[:, 0]]                          # [mm2] enabled cross sections with no entry in this database
    if unavailable:
        warnings.warn(f'MV connection: cable sections {unavailable} mm2 are not available as {material} '
                      f'{core} core cables at {V} kV, and have been excluded from the sizing')

    # keep only the rows whose cross section is both enabled and available
    mask = np.isin(full_database[:, 0], allowed)                                                        # [-] boolean mask over the database rows
    database = full_database[mask, :]                                                                   # [various] final database

    # without at least one usable cross section no sizing is possible
    if database.shape[0] == 0:
        raise ValueError('MV connection: no cable section available for the sizing. Please enable at least '
                         'one cross section compatible with the selected material, voltage and core')

    return database


def GetCableProp(Prop,
                 Material,
                 V,
                 Area,
                 Core,
                 Lay,
                 Duct,
                 ):
    """
    Returns one specific property of one specific cable, read from the database built by GetCableDatabase.

    Args:
        Prop (str):           [-] property required, one of:
                                  'R' - max conductor resistance @ 90 degC
                                  'C' - capacitance 'C' @ 50 Hz
                                  'L' - inductance 'L' @ 50 Hz, picked according to the laydown method
                                  'A' - ampacity, picked according to the laydown method and the duct flag
        Material (str):       [-] cable material, either 'cu'/'copper' or 'alu'/'aluminum'. Not case sensitive
        V (float):            [kV] connection line rated phase-to-phase voltage
        Area (float):         [mm2] cable cross section whose property is required
        Core (str):           [-] cable construction, either 'mono' or 'multi'. Not case sensitive
        Lay (str):            [-] cable laydown method, either 'trefoil' or 'flat'. Not case sensitive
        Duct (bool):          [-] True if the cables are laid inside a duct, False if directly buried

    Returns:
        value (float): the requested property, in the database units, i.e.:
                       Prop='R' -> [Ohm/km], Prop='C' -> [uF/km], Prop='L' -> [mH/km], Prop='A' -> [A]
    """

    # 1 - NORMALISE THE STRING ARGUMENTS
    prop = _assign_str_from_list(Prop, 'MV connection: cable property', ['r', 'c', 'l', 'a'])            # validate the requested property
    core = _assign_str_from_list(Core, 'MV connection: cable core', _VALID_CABLE_CORE)                   # validate the cable construction
    lay = _assign_str_from_list(Lay, 'MV connection: cable laydown', _VALID_CABLE_LAY)                   # validate the laydown method

    # 2 - LOAD THE DATABASE OF THE REQUESTED CROSS SECTION ONLY
    # the flags are all False but the one matching the requested cross section, so that the database returned
    # holds a single row
    flags = [sect == Area for sect in _VALID_MV_CABLE_SECT]                                             # [-] one flag per commercial cross section
    database = GetCableDatabase(Material, V, core, *flags)                                              # [various] single row database

    # 3 - PICK THE COLUMN HOLDING THE REQUESTED PROPERTY
    if prop == 'r':
        column = 1                                              # resistance, same column for both constructions
    elif prop == 'c':
        column = 2                                              # capacitance, same column for both constructions
    elif prop == 'l' and core == 'multi':
        column = 3                                              # multi core cables have one single inductance value
    elif prop == 'l' and lay == 'trefoil':
        column = 3                                              # single core cables in trefoil
    elif prop == 'l' and lay == 'flat':
        column = 4                                              # single core cables in flat
    elif prop == 'a' and core == 'multi':
        column = 5 if Duct else 4                               # multi core ampacity does not depend on the laydown method
    elif prop == 'a' and Duct:
        column = 7                                              # single core cables in duct are always laid in trefoil
    elif prop == 'a' and lay == 'trefoil':
        column = 6                                              # single core cables directly buried, in trefoil
    else:
        column = 5                                              # single core cables directly buried, in flat

    # 4 - RETURN THE VALUE
    return float(database[0, column])


# ---------------------------------------------------------------------------
#   DERATING FACTORS
# ---------------------------------------------------------------------------

def GetDerating(Nt,
                Area,
                Core,
                Lay,
                SoilR,
                GroundT,
                Depth,
                Spacing,
                Duct,
                ):
    """
    Computes the cumulative derating factor to be applied to the nominal cable ampacity, as the product of the
    five elementary derating factors (ground temperature, trench depth, soil resistivity, grouping of directly
    buried cables, grouping of cables in duct).

    Args:
        Nt (int):             [-] number of cable triads laid at the same depth
        Area (float):         [mm2] cable cross section
        Core (str):           [-] cable construction, either 'mono' or 'multi'. Not case sensitive
        Lay (str):            [-] cable laydown method, either 'trefoil' or 'flat'. Not case sensitive
        SoilR (float):        [K*m/W] soil thermal resistivity
        GroundT (float):      [degC] ground temperature
        Depth (float):        [m] trench depth
        Spacing (float):      [m] spacing between two adjacent cable triads
        Duct (bool):          [-] True if the cables are laid inside a duct, False if directly buried

    Returns:
        derating (float):     [-] cumulative derating factor, in the range [0, 1]
    """

    # the five factors are independent of each other, so the cumulative derating is simply their product
    return (getk1(GroundT)
            * getk2(Area, Depth, Core, Duct)
            * getk3(SoilR)
            * getk4(Lay, Spacing, Core, Duct, Nt)
            * getk5(Duct, Nt))


def getk1(GroundT):
    """
    Computes the ground temperature derating factor.

    Args:
        GroundT (float):      [degC] ground temperature

    Returns:
        k1 (float):           [-] ground temperature derating factor
    """

    # linear interpolation over the tabulated values. Outside the tabulated range np.interp holds the closest
    # tabulated value, which is a safe behaviour since the input check already bounds the temperature
    return float(np.interp(GroundT, _K1_TEMPERATURE, _K1_FACTOR))


def getk2(Area,
          Depth,
          Core,
          Duct,
          ):
    """
    Computes the trench depth derating factor.

    Args:
        Area (float):         [mm2] cable cross section
        Depth (float):        [m] trench depth
        Core (str):           [-] cable construction, either 'mono' or 'multi'. Not case sensitive
        Duct (bool):          [-] True if the cables are laid inside a duct, False if directly buried

    Returns:
        k2 (float):           [-] trench depth derating factor
    """

    # normalise the cable construction, so that the lookup key can be built safely
    core = _assign_str_from_list(Core, 'MV connection: cable core', _VALID_CABLE_CORE)

    # pick the installation case. Single core cables have two tables, split by cross section, because large
    # conductors dissipate more heat per metre of trench and are therefore more sensitive to the depth
    if core == 'multi':
        size = 'any'                                                            # multi core cables have one single table
    else:
        size = 'small' if Area <= _K2_AREA_THRESHOLD else 'large'               # single core cables are split at 185 mm2

    # build the lookup key and read the corresponding row of tabulated factors
    factors = _K2_FACTOR[('duct' if Duct else 'buried', core, size)]            # [-] tabulated factors for this installation case

    # linear interpolation over the trench depth
    return float(np.interp(Depth, _K2_DEPTH, factors))


def getk3(SoilR):
    """
    Computes the soil thermal resistivity derating factor.

    Args:
        SoilR (float):        [K*m/W] soil thermal resistivity

    Returns:
        k3 (float):           [-] soil thermal resistivity derating factor
    """

    # linear interpolation over the tabulated values. A resistivity below 1 K*m/W gives a factor above 1,
    # since the reference ampacity of the datasheet is given for a 1 K*m/W soil
    return float(np.interp(SoilR, _K3_RESISTIVITY, _K3_FACTOR))


def getk4(Lay,
          Spacing,
          Core,
          Duct,
          Nt,
          ):
    """
    Computes the grouping derating factor of directly buried cables, i.e. the mutual heating between the
    cable triads laid side by side at the same depth.

    Args:
        Lay (str):            [-] cable laydown method, either 'trefoil' or 'flat'. Not case sensitive
        Spacing (float):      [m] spacing between two adjacent cable triads
        Core (str):           [-] cable construction, either 'mono' or 'multi'. Not case sensitive
        Duct (bool):          [-] True if the cables are laid inside a duct, False if directly buried
        Nt (int):             [-] number of cable triads laid at the same depth

    Returns:
        k4 (float):           [-] grouping derating factor for directly buried cables
    """

    # cables laid in a duct are not derated here: their mutual heating is accounted for by k5 instead
    if Duct:
        return 1.0

    # normalise the string arguments, so that the lookup key can be built safely
    core = _assign_str_from_list(Core, 'MV connection: cable core', _VALID_CABLE_CORE)
    lay = _assign_str_from_list(Lay, 'MV connection: cable laydown', _VALID_CABLE_LAY)

    # the grouping tables stop at 6 triads. Beyond that, the derating is assumed to be the one of a full layer,
    # since 6 is the max number of triads that can be laid at the same depth
    n_triads = int(min(max(Nt, 1), _K4_MAX_TRIADS))                             # [-] number of triads, clipped to the tabulated range

    # read the row of tabulated factors of this triad count, then interpolate over the spacing
    factors = _K4_FACTOR[(lay, core)][n_triads - 1, :]                          # [-] tabulated factors for this installation case

    return float(np.interp(Spacing, _K4_SPACING, factors))


def getk5(Duct,
          Nt,
          ):
    """
    Computes the grouping derating factor of cables laid in a duct, i.e. the mutual heating between the cable
    triads laid side by side inside the same duct bank.

    Args:
        Duct (bool):          [-] True if the cables are laid inside a duct, False if directly buried
        Nt (int):             [-] number of cable triads laid at the same depth

    Returns:
        k5 (float):           [-] grouping derating factor for cables in duct
    """

    # directly buried cables are not derated here: their mutual heating is accounted for by k4 instead
    if not Duct:
        return 1.0

    # 6 triads, i.e. a full duct bank layer, have their own tabulated factor. More than 6 triads cannot be laid
    # at the same depth, so the derating of a full layer is held: the factor must not grow back with Nt
    if Nt >= _K4_MAX_TRIADS:
        return _K5_SIX_TRIADS

    # up to 5 triads, the in-duct grouping is read from the touching multi core trefoil table, shifted by one
    # position: cables in a duct are always in contact with the duct wall, so even a single triad is derated
    return getk4('trefoil', 0, 'multi', False, Nt + 1)



# ---------------------------------------------------------------------------
#   INPUT CHECK
# ---------------------------------------------------------------------------

def GetCableSizeOpt_input_check(Pdc,
                                Pac,
                                EqHours,
                                Lenght,
                                V,
                                Material,
                                MaxEloss=1.5,
                                Core='mono',
                                Lay='trefoil',
                                Cosphi=0.9,
                                SoilR=1.5,
                                GroundT=30,
                                Depth=0.8,
                                Spacing=0.15,
                                Duct=True,
                                MaxLoading=0.8,
                                BetaA=0.53,
                                BetaB=0.57,
                                f50=True,
                                f70=True,
                                f95=True,
                                f120=True,
                                f150=True,
                                f185=True,
                                f240=True,
                                f300=True,
                                f400=True,
                                f500=True,
                                f630=True,
                                iterations=3,
                                seed=None,
                                freq=50,
                                ):
    """
    Validates every argument of GetCableSizeOpt, before any calculation is started.

    Mandatory arguments (plant data) have no fallback value: an invalid or missing value raises an exception.
    Optional arguments (installation data) fall back on their standard value and raise a warning instead,
    so that a partially filled input set still produces a usable sizing.

    Args:
        Pdc (float):          [MWp] DC power. Used only for PV plants. For BESS, DC power equals AC power
        Pac (float):          [MVA] AC power at the injection point
        EqHours (float):      [kWh/kWp or h] plant yield (PV) or equivalent full load hours (BESS)
        Lenght (float):       [km] connection line length
        V (float):            [kV] connection line rated phase-to-phase voltage
        Material (str):       [-] cable material, either 'cu'/'copper' or 'alu'/'aluminum'. Not case sensitive
        MaxEloss (float):     [%] max admissible energy loss along the connection line
        Core (str):           [-] cable construction, either 'mono' or 'multi'. Not case sensitive
        Lay (str):            [-] cable laydown method, either 'trefoil' or 'flat'. Not case sensitive
        Cosphi (float):       [-] power factor along the connection line
        SoilR (float):        [K*m/W] soil thermal resistivity
        GroundT (float):      [degC] ground temperature
        Depth (float):        [m] trench depth
        Spacing (float):      [m] spacing between two adjacent cable triads
        Duct (bool):          [-] True if the cables are laid inside a duct, False if directly buried
        MaxLoading (float):   [-] max loading allowed on the cable. It is needed because the derating factors
                                  do not account for the short circuit current, so an extra safety margin is applied
        BetaA (float):        [-] alpha parameter of the Beta distribution of the plant power output
        BetaB (float):        [-] beta parameter of the Beta distribution of the plant power output
        f50 (bool):           [-] flag: allow the use of a 50 mm2 cable
        f70 (bool):           [-] flag: allow the use of a 70 mm2 cable
        f95 (bool):           [-] flag: allow the use of a 95 mm2 cable
        f120 (bool):          [-] flag: allow the use of a 120 mm2 cable
        f150 (bool):          [-] flag: allow the use of a 150 mm2 cable
        f185 (bool):          [-] flag: allow the use of a 185 mm2 cable
        f240 (bool):          [-] flag: allow the use of a 240 mm2 cable
        f300 (bool):          [-] flag: allow the use of a 300 mm2 cable
        f400 (bool):          [-] flag: allow the use of a 400 mm2 cable. Single core cables only
        f500 (bool):          [-] flag: allow the use of a 500 mm2 cable. Single core cables only
        f630 (bool):          [-] flag: allow the use of a 630 mm2 cable. Single core cables only
        iterations (int):     [-] number of Monte Carlo runs used to smooth the energy loss estimate
        seed (int or None):   [-] seed of the random number generator. Set an integer for reproducible results
        freq (float):         [Hz] grid frequency, used to compute the cable inductive reactance

    Returns:
        checked (tuple): the validated arguments, in the following order:
            Pdc [MWp], Pac [MVA], EqHours [h], Lenght [km], V [kV], Material [-], MaxEloss [%], Core [-],
            Lay [-], Cosphi [-], SoilR [K*m/W], GroundT [degC], Depth [m], Spacing [m], Duct [-],
            MaxLoading [-], BetaA [-], BetaB [-], flags [-] (np.ndarray of 11 bool, ordered as
            _VALID_MV_CABLE_SECT), iterations [-], seed [-], freq [Hz]
    """

    # PLANT DATA - mandatory arguments, no fallback value is allowed
    Pdc = _assign_numeric(Pdc, 'MV connection: Pdc', lo=0, unit='MWp')                                  # DC power must be a positive number
    Pac = _assign_numeric(Pac, 'MV connection: Pac', lo=0, unit='MVA')                                  # AC power must be a positive number
    EqHours = _assign_numeric(EqHours, 'MV connection: EqHours', lo=0, hi=8760, unit='h')               # equivalent hours cannot exceed the hours in one year
    Lenght = _assign_numeric(Lenght, 'MV connection: line length', lo=0, unit='km')                     # line length must be a positive number
    V = _assign_numeric(V, 'MV connection: line voltage', lo=0, hi=max(_VOLTAGE_CLASSES), unit='kV')    # voltage must fall inside the range covered by the cable database
    Material = _assign_str_from_list(Material, 'MV connection: cable material', _VALID_CABLE_MATERIAL)  # material must be one of the accepted aliases

    # normalise the material aliases, so that only 'cu' and 'alu' circulate inside the module
    Material = {'copper': 'cu', 'aluminum': 'alu'}.get(Material, Material)

    # aluminium cables are not manufactured at 36 kV, only copper ones are
    if Material == 'alu' and V > _VOLTAGE_CLASSES[-2]:
        raise ValueError(f'MV connection: aluminium cables are not available above {_VOLTAGE_CLASSES[-2]} kV. '
                         f'Please select copper cables or lower the connection voltage')

    # INSTALLATION DATA - optional arguments, an invalid value falls back on its standard value
    MaxEloss = _assign_numeric(MaxEloss, 'MV connection: MaxEloss', lo=0, hi=100, standard=1.5, unit='%')                                                # max admissible energy loss, as a percentage
    Core = _assign_str_from_list(Core, 'MV connection: cable core', _VALID_CABLE_CORE, standard='mono', custom_msg=' (single core cables)')              # cable construction
    Lay = _assign_str_from_list(Lay, 'MV connection: cable laydown', _VALID_CABLE_LAY, standard='trefoil', custom_msg=' (cables laid in trefoil)')       # cable laydown method
    Cosphi = _assign_numeric(Cosphi, 'MV connection: Cosphi', lo=0, hi=1, standard=0.9, unit='-')                                                        # power factor along the line
    SoilR = _assign_numeric(SoilR, 'MV connection: soil thermal resistivity', lo=_K3_RESISTIVITY[0], hi=_K3_RESISTIVITY[-1], standard=1.5, unit='K*m/W') # soil thermal resistivity, bounded by the k3 table
    GroundT = _assign_numeric(GroundT, 'MV connection: ground temperature', lo=_K1_TEMPERATURE[0], hi=_K1_TEMPERATURE[-1], standard=30, unit='degC')     # ground temperature, bounded by the k1 table
    Depth = _assign_numeric(Depth, 'MV connection: trench depth', lo=_K2_DEPTH[0], hi=_K2_DEPTH[-1], standard=0.8, unit='m')                             # trench depth, bounded by the k2 table
    Spacing = _assign_numeric(Spacing, 'MV connection: cable spacing', lo=_K4_SPACING[0], hi=_K4_SPACING[-1], standard=0.15, unit='m')                   # spacing between triads, bounded by the k4 table
    Duct = _assign_bool(Duct, 'MV connection: duct laying', standard=True, custom_msg=', thus cables are assumed to be laid in a duct')                  # duct or directly buried
    MaxLoading = _assign_numeric(MaxLoading, 'MV connection: MaxLoading', lo=0, hi=1, standard=0.8, unit='-')                                            # extra safety margin on the cable ampacity
    BetaA = _assign_numeric(BetaA, 'MV connection: BetaA', lo=0, hi=1, standard=0.53, unit='-')                                                          # alpha parameter of the Beta distribution
    BetaB = _assign_numeric(BetaB, 'MV connection: BetaB', lo=0, hi=1, standard=0.57, unit='-')                                                          # beta parameter of the Beta distribution

    # CABLE SIZE FLAGS - collected in a single array, ordered exactly as _VALID_MV_CABLE_SECT
    flags_in = [f50, f70, f95, f120, f150, f185, f240, f300, f400, f500, f630]
    flags = np.array([_assign_bool(flag, f'MV connection: f{sect} flag', standard=True, custom_msg=f', thus the {sect} mm2 cable is allowed')
                      for flag, sect in zip(flags_in, _VALID_MV_CABLE_SECT)], dtype=bool)

    # MONTE CARLO AND GRID SETTINGS
    iterations = _assign_numeric(iterations, 'MV connection: Monte Carlo iterations', lo=1, standard=3, unit='-', int_only=True)     # number of Monte Carlo runs
    seed = _assign_numeric(seed, 'MV connection: random seed', lo=0, standard=None, unit='-', int_only=True)                         # random seed, None means a non reproducible run
    freq = _assign_numeric(freq, 'MV connection: grid frequency', lo=0, standard=50, unit='Hz')                                      # grid frequency

    # pack everything in the documented order
    return (Pdc, Pac, EqHours, Lenght, V, Material, MaxEloss, Core, Lay, Cosphi, SoilR, GroundT,
            Depth, Spacing, Duct, MaxLoading, BetaA, BetaB, flags, iterations, seed, freq)


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