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
This module contains the mathematical tools needed for
reaction calculations.
"""

from numpy import (sin, cos, arange, append, array,
                   sqrt, power, interp, polyval)
from scipy.interpolate import make_interp_spline as spline
from numba import njit


def mesh(r_min, r_max, dr):
    """
    mesh: creates a mesh with odd number of mesh points
    """
    r_space = arange(r_min, r_max + dr, dr)
    if len(r_space) % 2 == 0:
        r_space = arange(r_min, r_max + 2*dr, dr)
    return r_space

@njit
def f_der1(r_mesh, f_mesh):
    """
    Calculates first order derivative of a one dimensional function
    using central difference method.
    :math: `df(r)/dr = (f(r + dr) - f(r - dr)) / dr / 2`
    """
    dr = r_mesh[1] - r_mesh[0]
    fr = lambda r: interp(r, r_mesh, f_mesh)
    dfr = (fr(r_mesh + dr) - fr(r_mesh - dr)) / dr / 2
    return dfr

@njit
def f_der2(r_mesh, f_mesh):
    """
    Calculates second order derivative of a one dimensional function
    using central difference method.
    :math: `d^2f(r)/dr^2 = (f(r + dr) - 2 f(r)  + f(r - dr)) / dr / dr`
    """
    dr = r_mesh[1] - r_mesh[0]
    fr = lambda r: interp(r, r_mesh, f_mesh)
    dfr = (fr(r_mesh + dr) + fr(r_mesh - dr) - 2 * fr(r_mesh)) / dr / dr
    return dfr


@njit
def j_hat_1(r):
    """
    Calculates j_hat_1 = j_n(1, r) * 3/r
    """
    r3 = power(r, 3)

    # keep away r = r_min ~ 0 (1e-18)
    # this may cause true_divide or ZeroDivision
    j_hat = 3*(sin(r) - r * cos(r))/r3
    j_hat[0] = 1.0
    return j_hat

@njit
def j_n(n, r):
    fr = r.copy()
    for i in range(*fr.shape,):
        fr[i] = j_n_scalar(n, r[i])
    return fr

@njit
def j_n_scalar(n, r):
    """
    Calculates j(n, r): spherical  bessel function
    this function might be improved or deleted later.
    """
    # keep away r = r_min ~ 0 (1e-18)
    # this may cause true_divide or ZeroDivision
    jsin = j_sin(n, r) * sin(r)
    jcos = j_cos(n, r) * cos(r) if n>0 else 0*r

    return jsin + jcos

@njit
def j_sin(n, r):
    """
    j(n, r): spherical  bessel function
           : = u(n, r) sin(r) + v(n, r) cos(r)
    Calculates u(n,r) function for n>=0
    """
    a, b = 0, 1/r
    for ni in range(1, n+1):
        a, b = b, (2*float(ni) - 1)*b/r - a

    return b

@njit
def j_cos(n, r):
    """
    j(n, r): spherical  bessel function
           : = u(n, r) sin(r) + v(n, r) cos(r)
    Calculates v(n,r) function for n>0
    """
    b = 1/r; a = b/r
    for ni in range(1, n + 2):
        a, b = b, (3 - 2*ni)*b/r - a

    return b * (-1)**( (n + 1) % 2)

def f_der1s(r_mesh, f_mesh, skip=None, f_mirror=-1):
    """
    Calculates first order derivative of a one dimensional function
    using central difference method.
    :math: `df(r)/dr = (f(r + dr) - f(r - dr)) / dr / 2`

    The function is smoothed using f_spline before and after the derivation.
    """
    f_smoothed = f_spline(r_mesh, f_mesh, skip=skip)
    f1 = f_spline(r_mesh, f_der1(r_mesh, f_smoothed), skip=skip, f_mirror=f_mirror)
    return f1


def f_der2s(r_mesh, f_mesh, skip=None, f_mirror=1):
    """
    Calculates first order derivative of a one dimensional function
    using central difference method.
    :math: `df(r)/dr = (f(r + dr) - f(r - dr)) / dr / 2`

    The function is smoothed using f_spline before and after the derivation.
    """
    f_smoothed = f_spline(r_mesh, f_mesh, skip=skip)
    f2 = f_spline(r_mesh, f_der2(r_mesh, f_smoothed), skip=skip, f_mirror=f_mirror)
    return f2


def f_spline(r_mesh, f_mesh, skip=None, f_mirror=1):
    """
    Creates new f_mesh using make_interp_spline (B-spline) from scipy.interpolate.
    This is useful to obtain derivatives of f_external by approximating to cubic spline.
    Otherwise, the derivatives are useless due to lack of precision.

    Smooth the function by skipping some part of the function.
    If skip = None then 2.5% will be skipped otherwise
    skip parameter with an integer value will define to skipping step.

    Reference:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.make_interp_spline.html
    """
    if skip == None:
        skip = int(len(r_mesh) * 2.5 / 100)

    # the function mirrored around r=0
    # therefore new range is between -r_mesh to r_mesh
    # f_mirror = 1 or -1 for even and odd functions, respectively.
    r_new = append(-r_mesh[::skip][1:][::-1], r_mesh[::skip])
    f_new = append(f_mirror * f_mesh[::skip][1:][::-1], f_mesh[::skip])

    # cubic spline (n = 3)
    f_smoothed = spline(r_new, f_new, k=3)
    return f_smoothed(r_mesh)


def f_extrapolate(r, r_min, r_max, f_poly_val):
    # ...
    # Engauge Digitizer
    # http://markummitchell.github.io/engauge-digitizer/
    r_poly = r[r>=r_min]
    r_poly = r_poly[r_poly<=r_max]
    f_poly = polyval(array(f_poly_val), r_poly)
    f_interp = interp(r, r_poly, f_poly)
    return f_spline(r, f_interp)

#..........................................................................#
#******** Tools for error analysis *********#
#..........................................................................#

def stdev(data):
    """
    data : 1D array
    mean : mean of the "data"
    stdev: standard deviation of the "data"

    returns:
    mean, stdev
    """
    n = len(data)
    mean = sum(data) / n
    distance = mean - data
    stdev = sqrt(sum(distance ** 2) / n)
    return mean, stdev

def mafe(yr, yc):
    """
    yr   : reference values (array)
    yc   : comapred values (array)
    n    : data number (integer)
    merr : mafe (Mean absolute fractional error) (float)
    serr : standard deviation of merr
    ferr : fractional individual errors

    ferr = abs(1 - yc/yr)
    merr = 1/n sum(ferr)

    returns:
    merr, serr, ferr

    Reference:
    https://en.wikipedia.org/wiki/Mean_absolute_percentage_error
    """
    ferr = abs(1 - yc / yr)
    merr, serr = stdev(ferr)
    return merr, serr, ferr

def mape(yr, yc):
    """
    yr    : reference values (array)
    yc    : comapred values (array)
    mperr : mape (Mean absolute percentage error) (float)
    sperr : standard deviation of mperr
    fperr : percentage individual errors


    mperr, sperr, fperr = 100 * mafe_err(yr, yc)

    returns:
    mperr, sperr, fperr

    Reference:
    https://en.wikipedia.org/wiki/Mean_absolute_percentage_error
    """
    return tuple([100 * m for m in mafe(yr, yc)])
    
def mse(yr, yc):
    """
    yr        : reference values (array)
    yc        : comapred values (array)
    mse_val   : mean squared error (mse)
    mse_dev   : standard deviation of mse
    distance2 : squared difference of yr and yc 
    
    returns:
    mse_val, mse_dev, distance2

    Reference:
    https://en.wikipedia.org/wiki/Mean_squared_error
    """
    distance = yr - yc
    distance2 = distance * distance
    mse_val, mse_dev = stdev(distance2)
    
    return mse_val, mse_dev, distance2

def rmse(yr, yc):
    """
    yr         : reference values (array)
    yc         : comapred values (array)
    rmse_val   : root mean squared error (rmse)
    rmse_dev   : standard deviation of rmse
    distance2  : squared difference of yr and yc 
    
    returns:
    rmse_val, rmse_dev, distance2

    Reference:
    https://en.wikipedia.org/wiki/Root-mean-square_deviation
    """
    
    mse_val, mse_dev, distance2 = mse(yr, yc)
    rmse_val = sqrt(mse_val)
    rmse_dev = sqrt(mse_dev)
    return rmse_val, rmse_dev, distance2

def wmse(yr, yc):
    """
    yr         : reference values (array)
    yc         : comapred values (array)
    wmse_val   : weighted mean squared error (wmse)
    wmse_dev   : standard deviation of wmse
    wdistance2 : weighted and squared difference of yr and yc 
    
    returns:
    wmse_val, wmse_dev, wdistance2

    Reference:
    Equation 38 in CPC 181 (1) (2010) 168–182.
    https://doi.org/10.1016/j.cpc.2009.09.007  
    """
    distance = yr - yc
    weight = yr + yc
    wdistance = distance/weight
    wdistance2 = wdistance * wdistance
    wmse_val, wmse_dev = stdev(wdistance2)
    
    return wmse_val, wmse_dev, wdistance2    