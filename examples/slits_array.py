"""Config file to simulate a membrane with a staggered lattice of rectangular slits"""

import numpy as np


# General parameters:
OUTPUT_FOLDER_NAME             = 'Array of slits'
NUMBER_OF_PHONONS              = 50
NUMBER_OF_TIMESTEPS            = 30000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 0.5e-12
T                              = 4.0
PLOTS_IN_TERMINAL              = False
OUTPUT_SCATTERING_MAP          = False
OUTPUT_RAW_THERMAL_MAP         = True
OUTPUT_TRAJECTORIES_OF_FIRST   = 30
NUMBER_OF_LENGTH_SEGMENTS      = 10


# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 100
NUMBER_OF_PIXELS_Y             = 100
NUMBER_OF_TIMEFRAMES           = 6


# Material parameters:
MEDIA                          = 'Si'
SPECIFIC_HEAT_CAPACITY         = 0.0176  # [J/kg/K] for Si at 4 K


# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 1200e-9
LENGTH                         = 2200e-9


# Phonon source:
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random_up")]


# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
HOLE_ROUGHNESS                 = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9


HOLES = []

# Staggered attice of holes:
size_x                         = 200e-9
size_y                         = 100e-9
period_x                       = 300e-9
period_y                       = 300e-9
first_hole_coordinate = 300e-9
number_of_periods_x = 5
number_of_periods_y = 3
for i in range(number_of_periods_y):
    for j in range(number_of_periods_x):
        x = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
        y = first_hole_coordinate + i * period_y * 2
        HOLES.append(RectangularHole(x, y, size_x, size_y))
for i in range(number_of_periods_y):
    for j in range(number_of_periods_x - 1):
        x = - 2 * period_x / 2 + j * period_x - period_x / 2
        y = first_hole_coordinate + period_y + i * period_y * 2
        HOLES.append(RectangularHole(x, y, size_x, size_y))
