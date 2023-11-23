"""
Modules provides scattering processes on various complex objects.
Each function determines whether the scattering should happen and
call corresponding function for scattering on corresponding primitive.
"""

from math import pi, cos, sin, tan, exp, sqrt, atan, asin, acos, ceil
from random import random
from numpy import sign, array, linspace, abs, where, diff
from numpy.linalg import norm
from scipy.optimize import bisect

from freepaths.config import cf
from freepaths.move import move
from freepaths.scattering_types import Scattering
from freepaths.scattering_primitives import *
from freepaths.scattering_parabolic import *
from freepaths.scatterers import *


def internal_scattering(ph, flight, scattering_types):
    """Check if the time passed since previous diffuse scattering event reached
    the time until an internal scattering event, and if yes, scatters randomly"""
    if flight.time_since_previous_scattering >= ph.time_of_internal_scattering:
        scattering_types.internal = random_scattering(ph)


def reinitialization(ph, scattering_types):
    """Re-thermalize (diffusely) phonon when it comes back to one of the hot sides"""
    x, y, _ = move(ph, cf.timestep)

    if cf.hot_side_position_bottom and y < 0:
        scattering_types.hot_side = horizontal_surface_up_scattering(ph, cf.side_wall_roughness, is_diffuse=True)

    if cf.hot_side_position_top and y > cf.length:
        scattering_types.hot_side = horizontal_surface_down_scattering(ph, cf.side_wall_roughness, is_diffuse=True)

    if cf.hot_side_position_right and x > cf.width/2:
        scattering_types.hot_side = vertical_surface_left_scattering(ph, cf.side_wall_roughness, is_diffuse=True)

    if cf.hot_side_position_left and x < - cf.width/2:
        scattering_types.hot_side = vertical_surface_right_scattering(ph, cf.side_wall_roughness, is_diffuse=True)


def scattering_on_circular_holes(ph, x0, y0, radius, scattering_types, x, y, z):
    """Check if a phonon strikes a circular hole and calculate the new direction"""

    # If phonon is inside the circle with given radius:
    if (x - x0)**2 + (y - y0)**2 <= radius**2:
        if y == y0: y += 1e-9 # Prevent division by zero
        tangent_theta = atan((x - x0)/(y - y0))
        scattering_types.holes = circle_outer_scattering(ph, tangent_theta, y, y0, cf.hole_roughness)


def scattering_on_rectangular_holes(ph, x0, y0, Lx, Ly, scattering_types, x, y, z):
    """Check if the phonon strikes a rectangular hole and calculate new direction"""

    # If the phonon is inside the rectangle:
    if (abs(x - x0) <= Lx / 2) and (abs(y - y0) <= Ly / 2):

        # Coordinate y of the intersection with the hole side:
        y1 = (y0 - y) + cos(ph.theta)*(Lx/2 - abs(x0 - x))/abs(sin(ph.theta))

        # Scattering on the left wall:
        if abs(y1) <= Ly/2 and x < x0:
            scattering_types.holes = vertical_surface_left_scattering(ph, cf.hole_roughness)

        # Scattering on the right wall:
        elif abs(y1) <= Ly/2 and x > x0:
            scattering_types.holes = vertical_surface_right_scattering(ph, cf.hole_roughness)

        # Scattering on the top wall:
        elif y > y0:
            scattering_types.holes = horizontal_surface_up_scattering(ph, cf.side_wall_roughness)

        # Scattering on the bottom wall:
        else:
            scattering_types.holes = horizontal_surface_down_scattering(ph, cf.side_wall_roughness)


def scattering_on_circular_pillars(ph, pillar, scattering_types, x, y, z):
    """Check if a phonon strikes a circular pillar and calculate new direction"""

    # Cone radius at a given z coordinate:
    radius = pillar.diameter/2# - (z - cf.thickness / 2) / tan(pillar.wall_angle)
    distance_from_pillar_center = sqrt((x - pillar.x)**2 + (y - pillar.y)**2)
    distance_from_pillar_center_original = sqrt((ph.x - pillar.x)**2 + (ph.y - pillar.y)**2)
    step = 2 * ph.speed * cf.timestep

    # If phonon crosses the pillar boundary. Third condition is to exclude all other pillars:
    if (distance_from_pillar_center >= radius and z > cf.thickness / 2 and distance_from_pillar_center < radius + step):

        # Calculate angle to the surface and specular scattering probability:
        tangent_theta = atan((x - pillar.x)/(y - pillar.y))
        scattering_types.pillars = circle_inner_scattering(ph, tangent_theta, y, pillar.y, cf.pillar_roughness)


def scattering_on_triangle_down_holes(ph, x0, y0, Lx, Ly, scattering_types, x, y, z):
    """Check if the phonon strikes a reverse triangular hole and calculate new direction after the scattering"""

    # Angle of the triangle:
    beta = atan(0.5 * Lx / Ly)

    # If phonon is inside the triangle:
    if (Ly/2 - (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2):

        # Scattering on the top wall of the triangle:
        if (ph.y > y0 + Ly / 2) and (abs(ph.theta) > pi / 2):
            scattering_types.holes = horizontal_surface_up_scattering(ph, cf.hole_roughness)

        # Scattering on the sidewalls of the triangle:
        else:
            scattering_types.holes = inclined_surfaces_down_scattering(ph, beta, x, x0, cf.hole_roughness)


def scattering_on_triangle_up_holes(ph, x0, y0, Lx, Ly, scattering_types, x, y, z):
    """Check if the phonon strikes a reverse triangular hole and calculate new direction"""

    # Angle of the triangle:
    beta = atan(0.5 * Lx / Ly)

    # If phonon is inside the triangle:
    if (Ly/2 + (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2):

        # Scattering on the bottom wall of the triangle:
        if (ph.y < y0 - Ly / 2) and (abs(ph.theta) < pi / 2):
            scattering_types.holes = horizontal_surface_down_scattering(ph, cf.hole_roughness)

        # Scattering on the sidewalls of the triangle:
        else:
            scattering_types.holes = inclined_surfaces_up_scattering(ph, beta, x, x0, cf.hole_roughness)


def scattering_on_half_triangle_up_holes(ph, x0, y0, Lx, Ly, is_right_half, scattering_types, x, y, z):
    """Check if the phonon strikes a reverse triangular hole and calculate new direction"""

    # Angle of the triangle:
    beta = atan(0.5 * Lx / Ly)

    # If phonon is inside the right side of the triangle:
    if is_right_half and (Ly/2 + (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2) and (x > x0):

            # Scattering on the bottom wall of the triangle:
            if (ph.y < y0 - Ly / 2) and (abs(ph.theta) < pi / 2):
                scattering_types.holes = horizontal_surface_down_scattering(ph, cf.hole_roughness)

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x < x0:
                scattering_types.holes = vertical_surface_left_scattering(ph, cf.hole_roughness)

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_up_scattering(ph, beta, x, x0, cf.hole_roughness)

    # If phonon is inside the left side of the triangle:
    if not is_right_half and (Ly/2 + (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2) and (x < x0):

            # Scattering on the bottom wall of the triangle:
            if (ph.y < y0 - Ly / 2) and (abs(ph.theta) < pi / 2):
                scattering_types.holes = horizontal_surface_down_scattering(ph, cf.hole_roughness)

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x > x0:
                scattering_types.holes = vertical_surface_right_scattering(ph, cf.hole_roughness)

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_up_scattering(ph, beta, x, x0, cf.hole_roughness)


def scattering_on_half_triangle_down_holes(ph, x0, y0, Lx, Ly, is_right_half, scattering_types, x, y, z):
    """Check if the phonon strikes a reverse triangular hole and calculate new direction"""

    # Angle of the triangle:
    beta = atan(0.5 * Lx / Ly)

    # If phonon is inside the right side of the triangle:
    if is_right_half and (Ly/2 - (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2) and (x > x0):

            # Scattering on the top wall of the triangle:
            if (ph.y > y0 + Ly / 2) and (abs(ph.theta) > pi / 2):
                scattering_types.holes = horizontal_surface_up_scattering(ph, cf.hole_roughness)

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x < x0:
                scattering_types.holes = vertical_surface_left_scattering(ph, cf.hole_roughness)

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_down_scattering(ph, beta, x, x0, cf.hole_roughness)

    # If phonon is inside the left side of the triangle:
    if not is_right_half and (Ly/2 - (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2) and (x < x0):

            # Scattering on the top wall of the triangle:
            if (ph.y > y0 + Ly / 2) and (abs(ph.theta) > pi / 2):
                scattering_types.holes = horizontal_surface_up_scattering(ph, cf.hole_roughness)

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x > x0:
                scattering_types.holes = vertical_surface_right_scattering(ph, cf.hole_roughness)

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_down_scattering(ph, beta, x, x0, cf.hole_roughness)


def scattering_on_sinus_wave(ph, box, sin_function, tolerance, bounds, thickness, scattering_types, bottom_points, top_points, xp, yp, zp):
    
    closest_distance = -1
    eval_bounds = None

    # first fast selection
    if box[0] < xp and xp < box[1] and box[2] < yp and yp < box[3]:
        # possibly inside the wave

        # check if the point is within one of the circles
        leftmost_point = sin_function(bounds[0])
        distance = numpy.linalg.norm(numpy.array(leftmost_point) - numpy.array([xp, yp]))
        if distance < thickness/2:
            closest_distance = distance
            closest_point = leftmost_point
            eval_bounds = (bounds[0], bounds[0]+thickness/2)
        
        rightmost_point = sin_function(bounds[1])
        distance = numpy.linalg.norm(numpy.array(rightmost_point) - numpy.array([xp, yp]))
        if distance < thickness/2:
            closest_distance = distance
            closest_point = rightmost_point
            eval_bounds = (bounds[1]-thickness/2, bounds[1])

        # if the point is within one of the circles check if it is not closer to the function close to the circle
        if eval_bounds is not None:
            eval_points = list(numpy.linspace(eval_bounds[0], eval_bounds[1], int(numpy.ceil((thickness/2)/tolerance))))[1:]
            function_values = [sin_function(i)[1] for i in eval_points]
            distances = [numpy.linalg.norm(numpy.array([xp, yp]) - numpy.array([i, u])) for i, u in zip(eval_points, function_values)]
            distance = min(distances)
            if distance < closest_distance:
                closest_distance = distance
                closest_index = numpy.where(distances == distance)[0]
                closest_point = (eval_points[closest_index[0]], function_values[closest_index[0]])

        # if the point was not in the circles do a fast check if the point is close to the function
        else:
            # find top and bottom limit
            top_index = numpy.argmin(numpy.abs(top_points[0] - xp))
            top_y = top_points[1,top_index]
            bottom_index = numpy.argmin(numpy.abs(bottom_points[0] - xp))
            bottom_y = bottom_points[1,bottom_index]

            if bottom_y < yp and yp < top_y:
                # point inside the slit
                
                # calculate distance from line
                eval_points = list(numpy.linspace(bounds[0], bounds[1], int(numpy.ceil((bounds[1]-bounds[0])/tolerance))))
                function_values = [sin_function(i)[1] for i in eval_points]
                distances = [numpy.linalg.norm(numpy.array([xp, yp]) - numpy.array([i, u])) for i, u in zip(eval_points, function_values)]
                distance = min(distances)
                if distance < thickness/2:
                    closest_distance = distance
                    closest_index = numpy.where(distances == distance)[0]
                    closest_point = (eval_points[closest_index[0]], function_values[closest_index[0]])

    if closest_distance != -1:
        # only scatter if is moving towars structure
        direction = (xp - ph.x, yp - ph.y)
        dot_product = numpy.dot(direction, (closest_point[0] - ph.x, closest_point[1] - ph.y))

        if 0 < dot_product:
            # should be same mechanics as scattering on circle
            if yp == closest_point[1]: y += 1e-9 # Prevent division by zero
            tangent_theta = atan((xp - closest_point[0])/(yp - closest_point[1]))
            scattering_types.holes = circle_outer_scattering(ph, tangent_theta, yp, closest_point[1], cf.hole_roughness)


def scattering_on_right_sidewall(ph, scattering_types, x, y, z):
    """Scatter phonon if it reached right side wall"""
    if x > cf.width/2:
        scattering_types.walls = vertical_surface_left_scattering(ph, cf.side_wall_roughness)


def scattering_on_left_sidewall(ph, scattering_types, x, y, z):
    """Scatter phonon if it reached left side wall"""
    if x < -cf.width/2:
        scattering_types.walls = vertical_surface_right_scattering(ph, cf.side_wall_roughness)


def scattering_on_top_sidewall(ph, scattering_types, x, y, z):
    """Check if the phonon hits top side wall and output new vector"""
    if y > cf.length:
        scattering_types.walls = horizontal_surface_down_scattering(ph, cf.side_wall_roughness)


def scattering_on_bottom_sidewall(ph, scattering_types, x, y, z):
    """Check if the phonon hits bottom side wall and output new vector"""
    if y < 0.0:
        scattering_types.walls = horizontal_surface_up_scattering(ph, cf.side_wall_roughness)


def floor_scattering(ph, scattering_types, x, y, z):
    """Check if the phonon hits the floor surface and calculate new angles"""
    if z < -cf.thickness/2:
        scattering_types.top_bottom = in_plane_surface_scattering(ph, cf.top_roughness)


def ceiling_scattering(ph, scattering_types, x, y, z):
    """Check if the phonon hits the ceiling surface and if this place has a pillar and output new vector"""
    if z > cf.thickness / 2:

        # Regular scattering if there are no pillars:
        if not cf.pillars:
            scattering_types.top_bottom = in_plane_surface_scattering(ph, cf.top_roughness)

        # If it is under a pillar, then scatter from the pillar's ceiling:
        else:
            for pillar in cf.pillars:
                distance_from_pillar_center = sqrt((x - pillar.x)**2 + (y - pillar.y)**2)
                is_under_pillar = (distance_from_pillar_center < pillar.diameter/2)
                if is_under_pillar and z > pillar.height + cf.thickness / 2:
                    scattering_types.top_bottom = in_plane_surface_scattering(ph, cf.top_roughness)
                    return
                elif is_under_pillar and z <= pillar.height + cf.thickness / 2:
                    return
                else:
                    pass

            # Regular scattering if phonon is not under any of the pillars:
            if ph.z < cf.thickness/2:
                scattering_types.top_bottom = in_plane_surface_scattering(ph, cf.top_roughness)


def surface_scattering(ph, scattering_types):
    """Check for a surface scattering on this step"""

    # Preliminary move to see if phonon would cross something:
    x, y, z = move(ph, cf.timestep)

    # Scattering on top and bottom surfaces:
    ceiling_scattering(ph, scattering_types, x, y, z)
    floor_scattering(ph, scattering_types, x, y, z)

    # Scattering on sidewalls:
    if cf.include_right_sidewall:
        scattering_on_right_sidewall(ph, scattering_types, x, y, z)
    if cf.include_left_sidewall:
        scattering_on_left_sidewall(ph, scattering_types, x, y, z)
    if cf.include_top_sidewall:
        scattering_on_top_sidewall(ph, scattering_types, x, y, z)
    if cf.include_bottom_sidewall:
        scattering_on_bottom_sidewall(ph, scattering_types, x, y, z)


    # Scattering on holes:
    if cf.holes:
        # Check for each hole and each hole type:
        for hole in cf.holes:
            if isinstance(hole, CircularHole):
                scattering_on_circular_holes(ph, hole.x, hole.y, hole.diameter/2, scattering_types, x, y, z)

            elif isinstance(hole, RectangularHole):
                scattering_on_rectangular_holes(ph, hole.x, hole.y, hole.size_x, hole.size_y, scattering_types, x, y, z)

            elif isinstance(hole, TriangularUpHole):
                scattering_on_triangle_up_holes(ph, hole.x, hole.y, hole.size_x, hole.size_y, scattering_types, x, y, z)

            elif isinstance(hole, TriangularDownHole):
                scattering_on_triangle_down_holes(ph, hole.x, hole.y, hole.size_x, hole.size_y, scattering_types, x, y, z)

            elif isinstance(hole, TriangularUpHalfHole):
                scattering_on_half_triangle_up_holes(ph, hole.x, hole.y, hole.size_x, hole.size_y, hole.is_right_half, scattering_types, x, y, z)

            elif isinstance(hole, TriangularDownHalfHole):
                scattering_on_half_triangle_down_holes(ph, hole.x, hole.y, hole.size_x, hole.size_y, hole.is_right_half, scattering_types, x, y, z)

            elif isinstance(hole, ParabolaTop):
                top_parabola_scattering(ph, hole, cf.side_wall_roughness, scattering_types, x, y, z)

            elif isinstance(hole, ParabolaBottom):
                bottom_parabola_scattering(ph, hole, cf.side_wall_roughness, scattering_types, x, y, z)

            elif isinstance(hole, SinusWave):
                scattering_on_sinus_wave(ph, hole.box, hole.sin_function, hole.tolerance, hole.bounds, hole.thickness, scattering_types, hole.bottom_points, hole.top_points, x, y, z)

            else:
                pass

            # If there was any scattering, then no need to check rest of the holes:
            if scattering_types.holes is not None:
                break

    # Check for each pillar:
    if cf.pillars:
        for pillar in cf.pillars:
            if isinstance(pillar, CircularPillar):
                scattering_on_circular_pillars(ph, pillar, scattering_types, x, y, z)

            # If there was any scattering, then no need to check other pillars:
            if scattering_types.pillars is not None:
                break

    # Correct angle if it became more than 180 degrees:
    ph.correct_angle()
