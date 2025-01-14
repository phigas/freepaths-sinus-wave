"""This module provides phonon class which generates and moves a phonon"""

from math import pi, asin, exp, log
from random import random, choice, randint
from numpy import sign
from scipy.constants import k, hbar
import numpy as np
import enum

from freepaths.config import cf
from freepaths.options import Distributions, Materials
import freepaths.move


class Phonon:
    """A phonon particle with various physical properties"""

    def __init__(self, material, branch_number=None, phonon_number=None):
        """Initialize a phonon by assigning initial properties"""
        self.branch_number = branch_number
        self.phonon_number = phonon_number
        self.x = None
        self.y = None
        self.z = None
        self.f = None
        self.phi = None
        self.theta = None
        self.speed = None
        self.first_timestep = randint(0, cf.number_of_virtual_timesteps)

        # Assigning initial properties of the phonon:
        if self.branch_number is None:
            self.branch_number = choice(range(3))
        self.f_max = max(material.dispersion[:, self.branch_number + 1])

        # Assigning initial coordinates and angles:
        source = choice(cf.phonon_sources)
        self.x, self.y, self.z = source.generate_coordinates()
        self.theta, self.phi = source.generate_angles()
        if cf.is_two_dimensional_material:
            self.phi = 0.0
            self.z = 0.0

        # Frequency is assigned based on Planckian distribution in case MFP sampling mode:
        if phonon_number is None:
            self.assign_frequency(material)

        # Otherwise, frequency is just asigned depending on the phonon number:
        else:
            f_upper = abs(material.dispersion[phonon_number + 1, self.branch_number + 1])
            f_lower = abs(material.dispersion[phonon_number, self.branch_number +1])
            self.f = (f_upper + f_lower) / 2

        self.assign_speed(material)
        self.assign_internal_scattering_time(material)

    @property
    def wavelength(self):
        """Calculate wavelength of the phonon"""
        return self.speed/self.f

    @property
    def is_in_system(self):
        """Checks if the phonon at this timestep did not reach the cold side.
        Depending on where user set cold sides, we check if phonon crossed that line"""
        is_inside_top = self.y < cf.length
        is_inside_bottom = self.y > 0
        is_inside_right = self.x < cf.width / 2.0
        is_inside_left = self.x > - cf.width / 2.0
        return ((not cf.cold_side_position_top or is_inside_top) and
                (not cf.cold_side_position_bottom or is_inside_bottom) and
                (not cf.cold_side_position_right or is_inside_right) and
                (not cf.cold_side_position_left or is_inside_left))

    def assign_frequency(self, material):
        """Assigning frequency with probability according to Planckian distribution"""

        # Frequency of the peak of the Plank distribution:
        f_peak = material.default_speed/(2*pi*hbar*material.default_speed/(2.82*k*cf.temp))

        # Density of states for f_max in Debye apparoximation:
        dos_max = 3*((2*pi*f_peak)**2)/(2*(pi**2)*(material.default_speed**3))

        # Bose-Einstein distribution for f_max:
        bose_einstein_max = 1/(exp((hbar*2*pi*f_peak)/(k*cf.temp)) - 1)

        # Peak of the distribution (needed for normalization)
        plank_distribution_max = dos_max*hbar*2*pi*f_peak*bose_einstein_max

        while True:                                                 # Until we obtain the frequency
            self.f = f_peak*5*random()                              # Let's draw a random frequency in the 0 - 5*f_peak range
            dos = 3*((2*pi*self.f)**2)/(2*(pi**2)*(material.default_speed**3))  # Calculate the DOS in Debye apparoximation
            bose_einstein = 1/(exp((hbar*2*pi*self.f)/(k*cf.temp))-1)           # And the Bose-Einstein distribution
            plank_distribution = dos*hbar*2*pi*self.f*bose_einstein             # Plank distribution

            # We take the normalized distribution at this frequency and draw a random number,
            # If the random number is lower than the distribution, we accept this frequency
            if random() < plank_distribution/plank_distribution_max and self.f < self.f_max:
                break

    def assign_speed(self, material):
        """Calculate group velocity dw/dk according to the frequency and polarization"""
        point_num = abs((np.abs(material.dispersion[:, self.branch_number + 1] - self.f)).argmin() - 1)
        d_w = 2*pi*abs(material.dispersion[point_num + 1, self.branch_number + 1]
                       - material.dispersion[point_num, self.branch_number + 1])
        d_k = abs(material.dispersion[point_num + 1, 0] - material.dispersion[point_num, 0])
        self.speed = d_w / d_k

    def assign_internal_scattering_time(self, material):
        """Determine relaxation time after which this phonon will undergo internal scattering"""

        if cf.use_gray_approximation_mfp:
            self.time_of_internal_scattering = cf.gray_approximation_mfp / self.speed
        else:
            omega = 2*pi*self.f

            # Depending on the material we assign different relaxation times:
            if material.name == Materials.Si:
                deb_temp = 152.0
                tau_impurity = 1/(2.95e-45 * (omega ** 4))
                tau_umklapp = 1/(0.95e-19 * (omega ** 2) * cf.temp * exp(-deb_temp / cf.temp))
                tau_internal = 1/((1/tau_impurity) + (1/tau_umklapp))

            elif material.name == Materials.SiC:    # Ref. Joshi et al, JAP 88, 265 (2000)
                # In SiC we also take into account 4 phonon scattering.
                deb_temp = 1200
                tau_impurity = 1/(8.46e-45 * (omega ** 4))
                tau_umklapp = 1/(6.16e-20 * (omega ** 2) * cf.temp * exp(-deb_temp / cf.temp))
                tau_4p = 1/(6.9e-23 * (cf.temp ** 2) * (omega ** 2))
                tau_internal = 1/((1/tau_impurity) + (1/tau_umklapp) + (1/tau_4p))

            elif material.name == Materials.Graphite:    # Ref. PRB 87, 115421 (2013)
                # In graphene we assume that material is perfect and no impurity scattering is present.
                deb_temp = 1000.0
                tau_umklapp = 1/(3.18e-25 * (omega ** 2) * (cf.temp ** 3) * exp(-deb_temp / (3*cf.temp)))
                tau_internal = 1/(1/tau_umklapp)

            else:
                raise ValueError('Specified material does not exist in the database.')

            # Final relaxation time is determined with some randomization [PRB 94, 174303 (2016)]:
            self.time_of_internal_scattering = -log(random()) * tau_internal

    def move(self):
        """Move a phonon in one timestep and return new coordinates"""
        self.x, self.y, self.z = freepaths.move.move(self, cf.timestep)

    def correct_angle(self):
        """Check if angles are out of the [-pi:pi] range and return them back to this range"""
        if abs(self.theta) > pi:
            self.theta -= sign(self.theta)*2*pi
