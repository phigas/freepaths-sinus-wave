"""Example to showcase the sinus wave hole type"""

OUTPUT_FOLDER_NAME             = "Sinus wave example"
NUMBER_OF_PHONONS              = 200
NUMBER_OF_TIMESTEPS            = 60000
TIMESTEP                       = 1e-12
T                              = 4
THICKNESS                      = 150e-9
WIDTH                          = 500e-9
LENGTH                         = 500e-9
OUTPUT_TRAJECTORIES_OF_FIRST   = 100

PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random_up")]

HOLES                          = [SinusWave(x=-200e-9, y=200e-9, len=300e-9, gap=100e-9, deviation=25e-9, thickness=100e-9, tolerance=1e-12)]
