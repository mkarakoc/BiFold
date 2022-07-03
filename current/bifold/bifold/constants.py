#! /usr/bin/python3

# Copyright (C) 2022 Mesut Karakoç <mesutkarakoc@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.

"""
This module provides the necessary physical constants.
"""

from numpy import sqrt, pi, power

# Close to zero
zero = 1e-10

# Mathematical constants
sqrt_pi = sqrt(pi)
pi2 = 2*pi
pi4 = 4*pi
pi_sqr = pi * pi
pi2_inv = 1 / 2 / pi_sqr # 1/(2*pi^2)
pi2_32 = 2 * power(pi, 3/2)


# Physical constants

# inverse of fine structure constant
# https://www.physics.nist.gov/cgi-bin/cuu/Value?alphinv
alpha_inv = 137.035999084 #137.035999084(21)

# Planck constant times Speed of light
# CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?hbcmevf
hbc = 197.3269804 # MeV.fm

#vacuum electric permittivity
# https://physics.nist.gov/cgi-bin/cuu/Value?ep0
epsilon0 = 8.8541878128e-12 # F/m = C^2/N.m^2 = C^2/J.m

# electron volt-joule relationship
# https://physics.nist.gov/cgi-bin/cuu/Value?Revj
one_eV_J = 1.602176634e-19 # Joule
one_J_eV = 1/one_eV_J # eV
one_J_MeV = one_J_eV * 1e-6 # MeV
one_m_fm = 1e15 # fm
one_fm_m = 1/one_m_fm # m

# elementary charge
# https://www.physics.nist.gov/cgi-bin/cuu/Value?e
electron_ch = 1.602176634e-19 # Colulomb

# square of unit charge
# long path:
# epsilon0_MeV = epsilon0/(one_J_MeV * one_m_fm)
# e2_coul = electron_ch*electron_ch
# e2_MeVfm = e2_coul/epsilon0_MeV/pi4  # e^2 <--- e^2/(4*pi*epsilon_0)
#
# short path:
# with this constant electrostatic p.e. of two positive point charges
# can be written as e^2/r [MeV] is equal to e^2/r/(4*pi*epsilon_0) [Joule]
e2 = hbc / alpha_inv # MeV.fm

# atomic mass constant energy equivalent in MeV
# https://www.physics.nist.gov/cgi-bin/cuu/Value?muc2mev
mu_c2 = 931.49410242 # MeV


# Units:
units = {'r_zero': None,'pi': None, 'sqrt_pi':None, 'pi2':None, 'pi4':None, 'pi_sqr':None, 'pi2_inv':None, 'pi2_32':None,
         'alpha_inv':None,
         'hbc':'MeV.fm', 'epsilon0':'F/m', 'one_eV_J':'J', 'one_J_eV':'eV', 'one_J_MeV':'MeV',
         'one_m_fm':'fm', 'one_fm_m':'m', 'electron_ch': 'C', 'e2':'MeV.fm',
         'mu_c2':'MeV'}

#
# This part prints out all the constants when only this module runs.
#
if __name__=='__main__':
    print('\nThe physical and mathematical constants:')
    for name, value in globals().copy().items():
        if '__' not in name and '<' not in str(value) and name!='units':
            unit = units[name] if units[name] != None else ''
            print(f'{name:18s} = {value:>22} {unit}')
