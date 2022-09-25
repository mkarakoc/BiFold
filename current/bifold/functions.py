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
This module provides the necessary functions for structure and reactions.
"""
from .simpson.integrals import simpson
from .constants import *

from numpy import (exp, loadtxt, interp, arange, append, array,
                   power, ndarray)
from io import StringIO

from html.parser import HTMLParser
from urllib.request import Request, urlopen
from numba import njit

class MLStripper(HTMLParser):
    """ clean html
    https://stackoverflow.com/a/925630
    https://stackoverflow.com/a/12982689
    """
    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.reset()
        self.fed = []

    def handle_data(self, d):
        self.fed.append(d)

    def handle_entityref(self, name):
        self.fed.append('&%s;' % name)

    def handle_charref(self, name):
        self.fed.append('&#%s;' % name)

    def get_data(self):
        return ''.join(self.fed)

def _strip_once(value):
    """
    Internal tag stripping utility used by strip_tags.
    """
    s = MLStripper()
    s.feed(value)
    s.close()
    return s.get_data()

def volumes_int(r, f, norm=None, L=0, name='', **kwargs):
    rL2 = power(r, L + 2)
    integrand2 = f * rL2
    integrand4 = integrand2 * rL2
    vol2 = pi4 * simpson(integrand2, r)
    vol4 = pi4 * simpson(integrand4, r)
    msr = vol4/vol2

    if norm == None:
        renorm = 1.0
    elif L == 0:
        renorm = norm / vol2
    else:
        renorm = norm / vol2 / pi4

    vol2 *= renorm
    vol4 *= renorm
    return_dict = {'name':name, 'L':L, 'norm':norm, 'renorm':renorm, 'vol2':vol2, 'vol4':vol4, 'msr':msr}
    if kwargs:
        return_dict.update(kwargs)
    return return_dict


@njit
def f_0(r, f):
    # correction for r ~ 0  [r<1e-10]
    # obtain f[0] using linear extrapolation

    tangent = (f[2] - f[1]) / (r[2] - r[1])
    intercept = (r[2]*f[1] - r[1]*f[2]) / (r[2] - r[1])
    return tangent * r[0] + intercept

@njit
def f__1(r,f):
    # correction for r = r_max
    # obtain f[r_max] using linear extrapolation
    return f_0(r[::-1], f[::-1])


def volumes(func):
    def inner(*args, **kwargs):
        func_name = func.__name__
        r = list(args)[0]
        fr = func(*args, **kwargs)

        if func_name!= 'f_dirac_delta':
            fr = append(f_0(r, fr), fr[1:])

        kkeys = kwargs.keys()
        #norm = None if not 'norm' in (kkeys:=kwargs.keys()) else kwargs['norm']
        norm = None if not 'norm' in kkeys else kwargs['norm']
        L = 0 if not 'L' in kkeys else kwargs['L']

        if func_name == 'f_rho_dd' or func_name == 'f_rho_bd':
            c = None if not 'c' in kkeys else kwargs['c']
            alpha = None if not 'c' in kkeys else kwargs['alpha']
            beta = None if not 'beta' in kkeys else kwargs['beta']
            gamma = None if not 'gamma' in kkeys else kwargs['gamma']
            n = None if not 'n' in kkeys else kwargs['n']
            fv = volumes_int(r, fr, norm=norm, L=L, name=func_name,
                             c=c, alpha=alpha, beta=beta, gamma=gamma, n=n)
        elif func_name == 'f_dirac_delta':
            fv = {'name': func_name, 'L': L, 'norm': norm, 'renorm':1.0, 'vol2': float(fr[0]), 'vol4': 0., 'msr': 0.}
        elif 'k_fermi' in func_name:
            Cs = None if not 'Cs' in kkeys else kwargs['Cs']
            fv = volumes_int(r, fr, norm=norm, L=L, name=func_name, Cs=Cs)
        else:
            fv = volumes_int(r, fr, norm=norm, L=L, name=func_name)


        fr *= fv['renorm']
        return keep_info([fr, [fv]])
    return inner


class keep_info:
    def __init__(self, value):
        if isinstance(value, int) or isinstance(value, float) or isinstance(value, ndarray):
            self.value = value
        else:
            self.value, self.info = value

    def __call__(self):
        return self.value

    def __neg__(self):
        return keep_info((-self.value, self.info))

    def __add__(self, other):
        if isinstance(other, int) or isinstance(other, float)  or isinstance(other, ndarray):
            return keep_info((self.value + other, [*self.info]))
        else:
            return keep_info((self.value + other.value, [*self.info, *other.info]))

    def __radd__(self, other):
        if isinstance(other, int) or isinstance(other, float)  or isinstance(other, ndarray):
            return keep_info((other + self.value, [*self.info]))
        else:
            return keep_info((other.value + self.value, [*self.info, *other.info]))

    def __sub__(self, other):
        if isinstance(other, int) or isinstance(other, float)  or isinstance(other, ndarray):
            return keep_info((self.value - other, [*self.info]))
        else:
            return keep_info((self.value - other.value, [*self.info, *other.info]))

    def __rsub__(self, other):
        if isinstance(other, int) or isinstance(other, float)  or isinstance(other, ndarray):
            return keep_info((other - self.value, [*self.info]))
        else:
            return keep_info((other.value-self.value, [*self.info, *other.info]))

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, float) or isinstance(other, ndarray):
            return keep_info((self.value * other, [*self.info]))
        else:
            return keep_info((self.value * other.value, [*self.info, *other.info]))

    def __rmul__(self, other):
        if isinstance(other, int) or isinstance(other, float) or isinstance(other, ndarray):
            return keep_info((other * self.value , [*self.info]))
        else:
            return keep_info((other.value * self.value, [*self.info, *other.info]))

    def __truediv__(self, other):
        if isinstance(other, int) or isinstance(other, float)  or isinstance(other, ndarray):
            return keep_info((self.value / other, [*self.info]))
        else:
            return keep_info((self.value / other.value, [*self.info, *other.info]))

    def __rtruediv__(self, other):
        if isinstance(other, int) or isinstance(other, float)  or isinstance(other, ndarray):
            return keep_info((other/self.value, [*self.info]))
        else:
            return keep_info((other.value/self.value , [*self.info, *other.info]))

    def __pow__(self, other):
        if isinstance(other, int) or isinstance(other, float)  or isinstance(other, ndarray):
            return keep_info((self.value ** other, [*self.info]))
        else:
            return keep_info((self.value ** other.value, [*self.info, *other.info]))

    def __rpow__(self, other):
        if isinstance(other, int) or isinstance(other, float)  or isinstance(other, ndarray):
            return keep_info((other ** self.value, [*self.info]))
        else:
            return keep_info((other.value ** self.value, [*self.info, *other.info]))

    def __eq__(self, other):
        if str(self.value) == str(other.value) and str(self.info) == str(other.info):
            return True
        else:
            return False

    def copy(self):
        return keep_info((self.value, [*self.info]))
        
    def pinfo(self):
        """prints info"""
        for key, val in self.info[0].items():
            print(f'{key:6s} = {val}')
    
    
#------------------------------------------------------------
@volumes
def f_exp_decay(r, V0, a, **kwargs):
    """Calculates Exponential function
    :math:`v(r) = V_0 \exp(-a r)`.

    Parameters
    ----------
    r     : f(r)
    V0    : depth or strength
    a     : diffuseness

    Returns
    -------
    float
        ``f`` values
    """
    return V0 * exp(-a*r)

@volumes
def f_yukawa(r, V0, a, b, **kwargs):
    """Calculates the 3 parameter Yukawa function
    :math:`v(r) = \frac{V_0}{b} \frac{\exp{-a r}}{r}`.

    Parameters
    ----------
    r     : f(r)
    V0, b : depth or strength
    a     : diffuseness

    Returns
    -------
    float
        ``f`` values
    """
    return V0/b * exp(-a*r)/r


#------------------------------------------------------------
@volumes
def f_2prm_fermi(r, V0, R, a, **kwargs):
    """Calculates the 2 parameter Fermi (or Woods-Saxon) function
    :math:`v(r) = \frac{V_0}{(1 + \exp{\frac{r-R}{a}})}`.

    Parameters
    ----------
    r  : f(r)
    V0 : depth or strength
    R  : radius or turning point
    a  : diffuseness

    Returns
    -------
    float
        ``f`` values
    """
    return V0/(1 + exp((r-R)/a))

@volumes
def f_3prm_fermi(r, V0, w, R, a, **kwargs):
    """Calculates the 3 parameter Fermi (or Woods-Saxon) function
        :math:`v(r) = (1 + w*r^2) \frac{V_0}{(1 + \exp{\frac{r-R}{a}})}`.

    Parameters
    ----------
    r  : f(r)
    V0 : depth or strength
    R  : radius or turning point
    a  : diffuseness
    w  : mostly defines shape of interior region

    Returns
    -------
    float
        ``f`` values
    """
    return  (1 + w*r*r) * V0/(1 + exp((r-R)/a))

@volumes
def f_2prm_gaussian(r, V0, a, **kwargs):
    """Calculates the 2 parameter Gaussian function
        :math:`v(r) = V_0 \exp{(r/a)^2}`.

    Parameters
    ----------
    r  : f(r)
    V0 : depth or strength
    a  : diffuseness

    Returns
    -------
    float
        ``f`` values
    """
    ra = r/a
    return V0*exp(-ra*ra)

@volumes
def f_3prm_gaussian(r, V0, w, a, **kwargs):
    """Calculates the 3 parameter Gaussian function
        :math:`v(r) = V_0 (1 + w r^2) \exp{(r/a)^2}`.

    Parameters
    ----------
    r  : f(r)
    V0 : depth or strength
    a  : diffuseness
    w  : mostly defines shape of interior region

    Returns
    -------
    float
        ``f`` values
    """
    ra = r/a
    return (1 + w*r*r) * V0*exp(-ra*ra)

@volumes
def f_sog(r, Ris, Qis, RP, Ze=1, **kwargs):
    # ATOMIC DATA AND NUCLEAR DATA TABLES 36,495536 (1987)
    # H.DE VRIES, C.W.DE JAGER, and C.DE VRIES
    # Sum-of-Gaussians
    gamma = sqrt(2/3) * RP
    gamma2 = gamma * gamma
    gamma3 = gamma2 * gamma
    rho_sog = 0
    for  Ri, Qi in zip(Ris, Qis):
        Ai =  Qi/ (pi2_32 * gamma3 * (1 + 2* Ri * Ri/gamma2))
        r_m = (r - Ri) / gamma
        r_p = (r + Ri) / gamma
        rho_sog += Ai * (exp(-r_m * r_m) + exp(-r_p * r_p))
    return Ze * rho_sog

@volumes
def f_rho_dd(r, rho, **kwargs):
    """Calculates density dependence
    :math:`f(\rho) = \rho \exp(- beta \rho)`.
    """
    beta = kwargs['beta'] if kwargs else 1
    return rho() * exp(-beta * rho())

@volumes
def f_rho_bd(r, rho, **kwargs):
    """Calculates density dependence
    :math:`f(\rho) = \rho \times \rho ^ n`.
    """
    n = kwargs['n'] if kwargs else 1
    return rho() * power(rho(), n)

@volumes
def f_dirac_delta(r, V0, **kwargs):
    """
    This concept taken from DFPOT
    THIS IS NOT AN ACTUAL DIRCAL DELTA FUNCTION!
    """
    #fv = {'name':'f_dirac_delta', 'L':L, 'norm':norm, 'vol2':float(V0), 'vol4':0., 'msr':0.}
    #return keep_info([append(V0, 0 * r[1:]), [fv]])
    return append(V0, 0 * r[1:])


def f_external(r, external_data, norm=None, L=0, data_format=0, msg=False):
    """ Creates an interpolated function from the proper
    data using numpy interp(r, r_data, f_data) .

    if external_data is a file name then the data is read
    from a file. if file_name is a dict  then the data is
    taken from the dict.

    Parameters
    ----------
    external_data  : it could be a file name or a dictionary
    data_format: the structure of the data in the file or the dict.
                 the data could be sorted vertically or horizontally
                 while keeping same sorting order.
    0: r f(r) [default format]
    1: dr
       r0
       f(r)
    2: n
       dr
       r0
       f(r)


    Returns
    -------
    array
        ``f`` values
    """
    if isinstance(external_data, str):
        with open(external_data, 'r') as f:
            e_data = f.read()
    elif isinstance(external_data, dict):
        e_name, = *external_data,
        e_data = external_data[e_name]
        external_data = e_name
    else:
        print('improper data in \'f_external()\': use a file or a dict with the correct format')
        quit()

    try:
        file = array([])
        for e_i in e_data.strip().split('\n'):
            if not e_i.startswith('#'):
                file = append(file, loadtxt(StringIO(e_i)))

        if data_format == 0:
            data_format_msg = \
                '\thorizontal:\tr\tf(r) --> \n\n\t-- OR --\n\n\tvertical:\n\tr\tf(r)\n\t|\t|\n\tV\tV'
            f_interp = interp(r, file[::2], file[1::2])
        elif data_format == 1:
            data_format_msg = \
                '\thorizontal:\tdr\tr0\tf(r) --> \n\n\t-- OR --\n\n\tvertical:\n\tdr\n\tr0\n\tf(r)\n\t|\n\tV'
            dr, r0 = file[:2]
            fr_file = file[2:]
            r_file = arange(r0, r0 + len(fr_file)*dr, dr)
            f_interp = interp(r, r_file, fr_file)
            if int(dr)>0:
                raise ValueError
        elif data_format == 2:
            data_format_msg = \
                '\thorizontal:\tn\tdr\tr0\tf(r) --> \n\n\t-- OR --\n\n\tvertical:\n\tn\n\tdr\n\tr0\n\tf(r)\n\t|\n\tV'
            n, dr, r0 = file[:3]
            r_file = arange(r0, r0 + n*dr, dr)
            f_interp = interp(r, r_file, file[3:])
        elif data_format == 3:
            data_format_msg = \
                '\thorizontal:\tr\tfn(r)\tfp(r) --> \n\n\t-- OR --\n\n\tvertical:\n\tr\tfn(r)\tfp(r)\n\t|\t|\t\t|\n\tV\tV\t\tV'
            f_interp = interp(r, file[::3], file[1::3] + file[2::3])
        else:
            print(f'unrecognized data format for: {external_data}')
            quit()
    except ValueError:
        print(f'check data format for: {external_data}')
        print(f'correct format must be:\n\t# title (optional)\n{data_format_msg} \nfor data_format = {data_format}.')
        quit()

    if msg==True:
        print(f'\nbe sure the correct format used for {external_data}:')
        print(data_format_msg)

    fr = f_interp
    fv = volumes_int(r, fr, norm=norm, L=L, name='f_external')
    fr *= fv['renorm']
    return keep_info([fr, [fv]])


def f_internet(r, url, norm=None, L=0, data_format=3):
    """Creates an interpolated function from a data file given in the "url" website.
    See f_external() for the reading formats of the data file.
    """
    user_agent = 'Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.0.7) Gecko/2009021910 Firefox/3.0.7'
    headers    = {'User-Agent': user_agent, }
    request    = Request(url, None, headers)
    response   = urlopen(request)
    data_html  = response.read().decode('utf8')
    data_list  = [di for di in _strip_once(data_html).split('\n') if di != '']
    data_range = [i for i, di in enumerate(data_list) if di.startswith('-')]
    data_text  = '\n'.join([di if data_range[1]<i<data_range[2] else '#'+di for i, di in enumerate(data_list)])

    fr = f_external(r, {'f_ripl':data_text}, norm=norm, L=L, data_format=data_format).value
    fv = volumes_int(r, fr, norm=norm, L=L, name='f_internet')
    fr *= fv['renorm']
    return keep_info([fr, [fv]])

def f_ripl(r, Z=8, A=16, norm=None, L=0, data_format=3):
    """Creates an interpolated function from a data file given in the RIPL website.
    See f_external() for the reading formats of the data file.
    """
    url = f'https://www-nds.iaea.org/cgi-bin/ripl_masses_nmd.pl?Z={Z}&A={A}'
    fr = f_internet(r, url, norm=norm, L=L, data_format=data_format).value
    fv = volumes_int(r, fr, norm=norm, L=L, name='f_ripl')
    fr *= fv['renorm']
    return keep_info([fr, [fv]])

