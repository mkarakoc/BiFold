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
Simpson's integrations.
"""

from ..time_check import timer
from numpy import (sum)
from numba import njit
from ..matematik import *

###########################
# Simpson integration with Numba
##########################
@njit
#@cc.export('j_simpson', 'f8(f8[:], f8[:])')
def simpson(f, r):
    """
    Calculates integrate of f function over r
    using Simpson 1/3 rule.
    """
    dr = r[1] - r[0]
    f_even = sum(f[:-1][2::2])
    f_odd  = sum(f[1::2])
    return (f[0] + 2*f_even + 4*f_odd + f[-1]) * dr/3

@njit
def fourier_with_simpson(f, r, q, n=0):
    fq = q.copy()
    fr2 = f * r*r
    for i in range(*fq.shape):
        qr = q[i] * r
        fq[i] = simpson(fr2 * j_n(n, qr), r)
    return fq

# @timer
@njit
def u_ex_with_simpson(dGRs, k, vnn_ex, R, s, n=0):
    u0_ex = R.copy()
    vnn_ex_s2 = vnn_ex * s*s
    for iR in range(*u0_ex.shape):
        ks = k[iR]*s
        u0_ex[iR] = simpson(dGRs[iR, :] * vnn_ex_s2 * j_n(n, ks), s)
    return u0_ex

# @timer
@njit
def fqs_with_simpson(fr2, r, g, s, q):
    fqs = q.reshape(-1, 1) * s.reshape(1, -1)
    _, ns = fqs.shape
    for j in range(ns):
        fqs[:, j] = fourier_with_simpson(fr2 * j_hat_1(g * s[j]), r, q)
    return fqs

# @timer
@njit
def gRs_with_simpson(dFqs, R, s, q, n=0):
    grs = R.reshape(-1, 1) * s.reshape(1, -1)
    nr, ns = grs.shape
    q2 = q * q
    for i in range(nr):
        qR = q * R[i]
        for j in range(ns):
            grs[i, j] = simpson( dFqs.T[j, :] * q2 * j_n(n, qR), q)
    return grs