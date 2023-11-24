"""
Module that provides various types of scattering objects like hole of different shapes.
User is expected to use one of these objects in the config file.
"""

import numpy

class CircularHole:
    """Shape of a circular hole"""
    def __init__(self, x=0, y=0, diameter=100e-9):
        self.x = x
        self.y = y
        self.diameter = diameter


class RectangularHole:
    """Shape of a rectangular hole"""
    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9):
        self.x = x
        self.y = y
        self.size_x = size_x
        self.size_y = size_y


class TriangularUpHole:
    """Shape of a triangular hole facing up"""
    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9):
        self.x = x
        self.y = y
        self.size_x = size_x
        self.size_y = size_y


class TriangularDownHole:
    """Shape of a triangular hole facing down"""
    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9):
        self.x = x
        self.y = y
        self.size_x = size_x
        self.size_y = size_y


class TriangularDownHalfHole:
    """Shape of a half triangular hole facing down"""
    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9, is_right_half=True):
        self.x = x
        self.y = y
        self.size_x = size_x
        self.size_y = size_y
        self.is_right_half = is_right_half


class TriangularUpHalfHole:
    """Shape of a half triangular hole facing up"""
    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9, is_right_half=True):
        self.x = x
        self.y = y
        self.size_x = size_x
        self.size_y = size_y
        self.is_right_half = is_right_half


class CircularPillar:
    """Shape of a circular pillar with inclined wall"""
    def __init__(self, x=0, y=0, diameter=200e-9, height=300e-9, wall_angle=numpy.pi/2):
        self.x = x
        self.y = y
        self.diameter = diameter
        self.height = height
        self.wall_angle = wall_angle


class ParabolaTop:
    """Shape of a parabolic wall"""
    def __init__(self, tip=0, focus=0):
        self.tip = tip
        self.focus = focus


class ParabolaBottom:
    """Shape of a parabolic wall"""
    def __init__(self, tip=0, focus=0):
        self.tip = tip
        self.focus = focus
        
class SinusWave:
    """Shape of thick sinusodial wave"""
    def __init__(self, x0=0, y0=0, len=400e-9, gap=50e-9, deviation=25e-9, thickness=75e-9, tolerance=5e-10, inverted=False):
        
        self.tolerance = tolerance
        self.thickness = thickness
        self.bounds = (x0+gap/2+thickness/2, x0+len+gap/2-thickness/2)
        if not inverted:
            self.sin_function = lambda x: numpy.array([x, y0-(numpy.cos((x-x0)*2*numpy.pi/(len+gap))-1)/2*deviation])
            # define box for fast phonon selection (xmin, xmax, ymin, ymax)
            self.box = (x0+gap/2, x0+gap/2+len, self.sin_function(self.bounds[0])[1]-thickness/2, y0+deviation+thickness/2)
            derivative_fun = lambda x: numpy.array([x*0+1, 2*numpy.pi*numpy.sin((x-x0)*2*numpy.pi/(len+gap))/(len+gap)/2*deviation])
        else:
            self.sin_function = lambda x: numpy.array([x, y0+(numpy.cos((x-x0)*2*numpy.pi/(len+gap))-1)/2*deviation])
            # define box for fast phonon selection (xmin, xmax, ymin, ymax)
            self.box = (x0+gap/2, x0+gap/2+len, y0-deviation-thickness/2, self.sin_function(self.bounds[0])[1]+thickness/2)
            derivative_fun = lambda x: numpy.array([x*0+1, -2*numpy.pi*numpy.sin((x-x0)*2*numpy.pi/(len+gap))/(len+gap)/2*deviation])

        xs_to_evaluate = numpy.linspace(self.bounds[0], self.bounds[1], int(numpy.ceil((len+gap)/tolerance)))
        tan_vector = derivative_fun(xs_to_evaluate)

        orth_vector = tan_vector[[1,0]]
        orth_vector[0] = orth_vector[0]*-1
        orth_vector = orth_vector/numpy.linalg.norm(orth_vector, axis=0)*thickness/2

        function_value = self.sin_function(xs_to_evaluate)
        self.bottom_points = function_value - orth_vector
        self.top_points = function_value + orth_vector
