# Copyright (C) 2022 Mesut KarakoÃ§ <mesutkarakoc@gmail.com>
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
This module contains the mathematical tools for
Filon's integrations.
"""

from ..time_check import timer
from numpy import (sin, cos, sum)
from numba import njit
from ..matematik import *

@njit
def f_0(r, f):
    # correction for r ~ 0  [r<1e-10]
    # obtain f[0] using linear extrapolation

    tangent = (f[2] - f[1]) / (r[2] - r[1])
    return f[1] - tangent*r[1]

########################################################
########################################################
from numpy import mod, pi

@njit
def sinm(x):
    return sin(mod(x, 2 * pi))

@njit
def cosm(x):
    return cos(mod(x, 2 * pi))

@njit
def alpha_beta_gamma(theta):
    t = theta

    t2 = t * t
    t3 = t2 * t
    sinmt = sinm(t)
    cosmt = cosm(t)

    alpha = 1 / t + sinmt * cosmt / t2 - 2 * sinmt * sinmt / t3
    beta = 2 * ((1 + cosmt * cosmt) / t2 - 2 * sinmt * cosmt / t3)
    gamma = 4 * (sinmt / t3 - cosmt / t2)

    return alpha, beta, gamma

@njit
def s_even(f, x, t):
    return sum(f * sinm(t * x)) - (f[-1] * sinm(t * x[-1]) + f[0] * sinm(t * x[0])) / 2

@njit
def s_odd(f, x, t):
    return sum(f * sinm(t * x))

@njit
def filon(g, x, t):
    f = g * x / t

    f[0] = f_0(x, f)

    dx = x[1] - x[0]
    alpha, beta, gamma = alpha_beta_gamma(t * dx)

    int_alpha = alpha * (f[0] * cosm(t * x[0]) - f[-1] * cosm(t * x[-1]))
    int_beta = beta * s_even(f[::2], x[::2], t)
    int_gamma = gamma * s_odd(f[1:][::2], x[1:][::2], t)
    return (int_alpha + int_beta + int_gamma) * dx

@njit
def fourier_with_filon(f, r, q, n=0):
    fq = q.copy()
    for i in range(*fq.shape):
        fq[i] = filon(f, r, q[i])
    fq[0] = f_0(q, fq)
    return fq

# @timer
@njit
def u_ex_with_filon(dGRs, k, vnn_ex, R, s, n=0):
    u0_ex = R.copy()
    for iR in range(*u0_ex.shape):
        u0_ex[iR] = filon(dGRs[iR, :] * vnn_ex, s, k[iR])
    u0_ex[0] = f_0(R, u0_ex)
    return u0_ex


# @timer
@njit
def fqs_with_filon(fr2, r, g, s, q):
    fqs = q.reshape(-1, 1) * s.reshape(1, -1)
    _, ns = fqs.shape
    for j in range(ns):
        fqs[:, j] = fourier_with_filon(fr2 * j_hat_1(g * s[j]), r, q)
    return fqs

# @timer
@njit
def gRs_with_filon(dFqs, R, s, q, n=0):
    grs = R.reshape(-1, 1) * s.reshape(1, -1)
    nr, ns = grs.shape
    for i in range(nr):
        for j in range(ns):
            grs[i, j] = filon( dFqs.T[j, :], q, R[i])
    return grs
########################################################
########################################################
