"""Config file to simulate a simple nanowire at 4 K"""

import numpy as np


# General parameters:
OUTPUT_FOLDER_NAME             = 'Nanowire'
NUMBER_OF_PHONONS              = 2000
NUMBER_OF_TIMESTEPS            = 30000
T                              = 4.0


# Multiprocessing
NUMBER_OF_PROCESSES = 10


# Material parameters:
MEDIA                          = 'Si'
SPECIFIC_HEAT_CAPACITY         = 0.0176  # [J/kg/K] for Si at 4 K
#SPECIFIC_HEAT_CAPACITY        = 714     # [J/kg/K] for Si at 300 K


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 200e-9
LENGTH                         = 2000e-9


# Map & profiles parameters:
pixel_size = 30e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PHONONS          = False


# Simulation time parameters:
TIMESTEP                       = 1.0e-12
total_simulation_time = 500e-9 # This should be at least a couple times the initialization time
NUMBER_OF_VIRTUAL_TIMESTEPS    = int(total_simulation_time / TIMESTEP)
initialization_time = 75e-9 # This should be set so that it is bigger than most phonons travel times
INITIALIZATION_TIMESTEPS       = int(initialization_time / TIMESTEP)
NUMBER_OF_INITIALIZATION_TIMEFRAMES = 3


# Phonon source:
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random_up")]
