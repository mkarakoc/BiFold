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
This module contains the interactions.
"""

from .functions import *
from numpy import polyval, array

@volumes
def v_coulomb(r):
    """Coulomb potential of two point charges"""
    return e2/r

#..........................................................................#
#**************************** M3Y - Reid/Paris ****************************#
#..........................................................................#
def v_m3y_reid_d(r, L=0):
    # Physics Letters B 342(1995) 6-12
    # D.T. Khoa, W. von Oertzen.
    # M3Y-Reid: Eqs. 1a and 1b

    # Physics Reports 285 (1997) 143-243
    # M.E. Brandan, G.R. Satchler.

    # Nuclear Physics A284 (1977) 399 - 419
    # G.Bertsch, J.Borysowicz, H.McManus, and W.G.Love.
    # Physics Reports, 55(3), 183–254.
    # Satchler, G.R., & Love, W.G.(1979).
    return f_yukawa(r, +7999.00, 4.0, 4.0, L=L) + \
           f_yukawa(r, -2134.25, 2.5, 2.5, L=L)

def v_m3y_reid_ex_zr(r, e_lab, a_proj, L=0):
    return f_dirac_delta(r, v_j00_m3y_reid(e_lab, a_proj), L=L)

def v_m3y_reid_ex_fr(r, L=0):
    # Physics Letters B 342(1995) 6-12
    # D.T. Khoa, W. von Oertzen.
    # M3Y-Reid: Eqs. 1a and 1b

    # Physics Reports 285 (1997) 143-243
    # M.E. Brandan, G.R. Satchler.

    return f_yukawa(r, +4631.38, 4.0, 4.0) + \
           f_yukawa(r, -1787.13, 2.5, 2.5) + \
           f_yukawa(r, -7.8474, 0.7072, 0.7072)


def v_j00_m3y_reid(e_lab, a_proj):
    """

    THERE IS A DIFFERENT VERSION OF THIS PSEUDO POTENTIAL?
    THINK IF YOU WOULD LIKE TO IMPLEMENT THIS!

    # Physics Reports 285 (1997) 143-243
    # M.E. Brandan, G.R. Satchler.

    Love, W. G., & Owen, L. W. (1975).
    Exchange effects from realistic interactions in the reformulated optical model.
    Nuclear Physics A, 239(1), 74–82.
    doi:10.1016/0375-9474(75)91133-1

    Satchler, G. R., & Love, W. G. (1979).
    Folding model potentials from realistic interactions for heavy-ion scattering.
    Physics Reports, 55(3), 183–254.
    doi:10.1016/0370-1573(79)90081-4

    Khoa, D. T., Satchler, G. R., & von Oertzen, W. (1995).
    Folding analysis of the elasticLi6+12C scattering: Knock-on exchange
    effects, energy dependence, and dynamical polarization potential.
    Physical Review C, 51(4), 2069–2084.
    doi:10.1103/physrevc.51.2069
    """
    return -276 * (1 - 0.005*e_lab/a_proj) # MeV.fm^3

def v_m3y_paris_d(r, L=0):
    # Physics Letters B 342(1995) 6 -12
    # D.T. Khoa, W. von Oertzen.
    # M3Y-Paris: Eqs. 2a and 21b

    # Physics Reports 285 (1997) 143-243
    # M.E. Brandan, G.R. Satchler.

    # Physics Reports 285 (1997) 143-243
    # M.E. Brandan, G.R. Satchler.
    # Nuclear Physics A398 (1983) 269 - 278
    # N.ANANTARAMAN, H.TOKI and G.F.BERTSCH
    # Nuclear Physics A668 (2000) 3 - 41
    # D.T. Khoa, G.R. Satchler.
    return f_yukawa(r, +11061.625, 4.0, 4.0, L=L) + \
           f_yukawa(r,  -2537.500, 2.5, 2.5, L=L)

def v_m3y_paris_ex_zr(r, e_lab, a_proj, L=0):
    return f_dirac_delta(r, v_j00_m3y_paris(e_lab, a_proj), L=L)

def v_m3y_paris_ex_fr(r, L=0):
    # Physics Letters B 342(1995) 6 -12
    # D.T. Khoa, W. von Oertzen.
    # M3Y-Paris: Eqs. 2a and 21b

    # Physics Reports 285 (1997) 143-243
    # M.E. Brandan, G.R. Satchler.

    # Nuclear Physics A398 (1983) 269 - 278
    # N.ANANTARAMAN, H.TOKI and G.F.BERTSCH
    # Nuclear Physics A668 (2000) 3 - 41
    # D.T. Khoa, G.R. Satchler.
    return f_yukawa(r, -1524.25, 4.0, 4.0) + \
           f_yukawa(r,  -518.75, 2.5, 2.5) + \
           f_yukawa(r, -7.8474, 0.7072, 0.7072)

def v_j00_m3y_paris(e_lab, a_proj):
    """
    # Physics Reports 285 (1997) 143-243
    # M.E. Brandan, G.R. Satchler.

    Love, W. G., & Owen, L. W. (1975).
    Exchange effects from realistic interactions in the reformulated optical model.
    Nuclear Physics A, 239(1), 74–82.
    doi:10.1016/0375-9474(75)91133-1

    Satchler, G. R., & Love, W. G. (1979).
    Folding model potentials from realistic interactions for heavy-ion scattering.
    Physics Reports, 55(3), 183–254.
    doi:10.1016/0370-1573(79)90081-4

    Khoa, D. T., Satchler, G. R., & von Oertzen, W. (1995).
    Folding analysis of the elasticLi6+12C scattering: Knock-on exchange
    effects, energy dependence, and dynamical polarization potential.
    Physical Review C, 51(4), 2069–2084.
    doi:10.1103/physrevc.51.2069
    """
    return -590 * (1 - 0.002*e_lab/a_proj) # MeV.fm^3
#..........................................................................#
#**************************** M3Y - Reid/Paris ****************************#
#..........................................................................#

#..........................................................................#
#****************** Density Dependent M3Y - Reid [DDM3Y] ******************#
#..........................................................................#
def v_ddm3y_reid_cab(e_lab, a_proj):
    """
    Original DDM3Y is defined by Ref. [1].
    Calculates the original energy dependent DDM3Y factors: C(E), alpha(E), beta(E).
    They are obtained by polynomial fitting of the Fig 1. of the Ref. [2] using numpy poly.
    Fig 1. is digitized using Engauge Digitizer Version 12.1.

    References:
         1. A.M. KOBOS, et al, Nuclear Physic A384 (1982) 65-87.
         2. M. EL-AZAB FARID and G. R. SATCHLER, Nuclear Physics A438 (1985) 525-535.

    """

    # latest digitization
    poly_alpha = array([+1.36544e-14, -5.24526e-12, +8.09459e-10,
                        -6.31517e-08, +2.53649e-06, -4.56887e-05,
                        +0.000239806, -0.001084510, +0.075644400, +3.58688])

    # previous digitization
    poly_beta = array([+1.05558e-14, -4.74991e-12, +8.94520e-10,
                       -9.10470e-08, +5.38698e-06, -0.000184957,
                       +0.003462820, -0.029982300, -0.055198500,
                       +11.683100000])

    poly_c = array([+1.14377e-17, -2.23584e-14, +7.57500e-12,
                    -1.10286e-09, +8.11007e-08, -2.93878e-06,
                    +3.81273e-05, +0.000327325, -0.015440300,
                    +0.527615])

    ea = e_lab/a_proj
    if  not 2.5<ea<90.5:
        print(f'e_lab/a_proj ={ea:9.4f} MeV/nucleon must be between 2.5 MeV - 90.5 MeV per nucleon')
        quit()
    return polyval(poly_c, ea), polyval(poly_alpha, ea), polyval(poly_beta, ea)
#..........................................................................#
#****************** Density Dependent M3Y - Reid [DDM3Y] ******************#
#..........................................................................#

#..........................................................................#
#******** Density Dependent M3Y - Reid/Paris [B/C/D-DM3Y: Xdm3yn] *********#
#..........................................................................#
def v_xdm3yn_cabgn(dd_name='bdm3y1', vnn_name='reid'):
    """
    Xdm3yn: X stands for d, b, c and n is an integer
    X=d, n=1 --> ddm3y1
    X=b, n=1 --> bdm3y1
    X=c, n=1 --> bdm3y1 etc.
    F(rho) = C [1 + alpha exp(-beta rho) - gamma rho ^ n]
    K [MeV]: ...
    """
    dd_name = dd_name.lower()
    vnn_name = vnn_name.lower()
    params_dict = {
        'reid':{        #C      alpha     beta   gamma   n    K [MeV]  Implemented
            'dim3y':  (  None,    None,   None,   None, None, None,    True),
            'ddm3y1': (0.2845,  3.6391, 2.9605, 0.0000,  0.0,  171,    True),
            'bdm3y0': (1.3827,   0.000, 0.0000, 1.1135,  2/3,  232,    False),
            'bdm3y1': (1.2253,   0.000, 0.0000, 1.5124,  1.0,  232,    True),
            'bdm3y2': (1.0678,   0.000, 0.0000, 5.1069,  2.0,  354,    True),
            'bdm3y3': (1.0153,   0.000, 0.0000, 21.073,  3.0,  475,    True)
        },
        'paris': {      # C     alpha    beta     gamma   n    K [MeV]  Implemented
            'dim3y':  (  None,    None,   None,    None, None,  None,   True),
            'ddm3y1': (0.2963,  3.7231, 3.7384,  0.0000,  0.0,   176,   True),
            'bdm3y1': (1.2521,  0.0000,  0.000,  1.7452,  1.0,   270,   True),
            'bdm3y2': (1.0664,  0.0000,  0.000,  6.0296,  2.0,   418,   True),
            'bdm3y3': (1.0045,  0.0000,  0.000, 25.1150,  3.0,   566,   True),
                         # C    alpha   beta  gamma  n, K [MeV] Imlemented
            'cdm3y1': (0.3429, 3.0232, 3.5512, 0.5,  1,   188,    True),
            'cdm3y2': (0.3346, 3.0357, 3.0685, 1.0,  1,   204,    True),
            'cdm3y3': (0.2985, 3.4528, 2.6388, 1.5,  1,   217,    True),
            'cdm3y4': (0.3052, 3.2998, 2.3180, 2.0,  1,   228,    True),
            'cdm3y5': (0.2728, 3.7367, 1.8294, 3.0,  1,   241,    True),
            'cdm3y6': (0.2658, 3.8033, 1.4099, 4.0,  1,   252,    True)
        }}
    return params_dict[vnn_name][dd_name]
#..........................................................................#
#******** Density Dependent M3Y - Reid/Paris [B/C/D-DM3Y: Xdm3yn] *********#
#..........................................................................#
