"""Example to showcase the sinus wave hole type"""

# set simulation parameters
OUTPUT_FOLDER_NAME             = "Sinus wave example structure"
NUMBER_OF_PHONONS              = 100
NUMBER_OF_TIMESTEPS            = 60000
TIMESTEP                       = 1e-12
T                              = 50
OUTPUT_TRAJECTORIES_OF_FIRST   = 100

# set slit dimenstions
slit_thickness      = 75e-9
deviation           = 50e-9
slit_length         = 300e-9
x_gap               = 50e-9
y_gap               = x_gap

# define simulation size
padding             = 100e-9
x_num               = 2
y_num               = 4

# set tolerance for calculations
tolerance           = 1e-9

# calculate domain size
THICKNESS                      = 150e-9
WIDTH                          = x_num*(slit_length+x_gap) # in x direction
LENGTH                         = 2*padding+y_num*(slit_thickness+deviation+y_gap)

# add sources
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random_up")]

# build holes
HOLES = []
for i_y in range(y_num):
    x_pos = -WIDTH/2 - (slit_length+x_gap)/2
    y_pos = padding+(slit_thickness+deviation+y_gap)*i_y + deviation + y_gap + slit_thickness/2
    
    # HOLES.append(SinusWave(x0=x_pos, y0=y_pos, len=slit_length, gap=x_gap, deviation=deviation, thickness=slit_thickness, tolerance=tolerance, inverted=True))
    
    for i_x in range(x_num):
        x_pos = -WIDTH/2+(slit_length+x_gap)*i_x
        y_pos = padding+(slit_thickness+deviation+y_gap)*i_y
        
        HOLES.append(SinusWave(x0=x_pos, y0=y_pos, len=slit_length, gap=x_gap, deviation=deviation, thickness=slit_thickness, tolerance=tolerance, inverted=False))
        
        x_pos += (slit_length+x_gap)/2
        y_pos += deviation + y_gap + slit_thickness/2
        
        HOLES.append(SinusWave(x0=x_pos, y0=y_pos, len=slit_length, gap=x_gap, deviation=deviation, thickness=slit_thickness, tolerance=tolerance, inverted=True))
