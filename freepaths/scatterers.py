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
    def __init__(self, x0=0, y0=0, a=400e-9, gap=50e-9, deviation=25e-9, thickness=75e-9, tolerance=5e-10):
        # assert dir in ('x+', 'x-', 'y+', 'y-')
        
        self.tolerance = tolerance
        self.thickness = thickness
        self.sin_function = lambda x: numpy.array([x, y0-(numpy.cos((x-x0)*2*numpy.pi/a)-1)/2*deviation])
        self.bounds = (x0+gap/2+thickness/2, x0+a-gap/2-thickness/2)
        self.grad_function = lambda x, x_p, y_p: 2*(x-x_p)   +   (  2*numpy.pi*  numpy.sin(2*numpy.pi*(x-x0)/a)  *  (y0-y_p+( -numpy.cos((2*numpy.pi*(x-x))/a)+1 )/(2*deviation))  )/(  deviation*a  )
        
        # define box for fast phonon selection (xmin, xmax, ymin, ymax)
        self.box = (x0+gap/2, x0-gap/2+a, self.sin_function(self.bounds[0])[1]-thickness/2, y0+deviation+thickness/2)
        
        # not implemented
        # # coordinate transformer to handle different orientations
        # if dir == 'x+':
        #     self.transform = lambda x,y: (x, y)
        # elif dir == 'x-':
        #     self.transform = lambda x,y: (x, y-(y0-y))
        # elif dir == 'y+':
        #     self.transform = lambda x,y: (x, y))
        # elif dir == 'y-':
        #     self.transform = lambda x,y: ()
        
        
