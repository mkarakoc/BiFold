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
This module calculates the double folding potentials with Simpson's integration.
"""

from ..graph_tools import *
from ..print_tools import *
from ..constants import *
from ..functions import *
from .matematik import *
from ..interactions import *

# Coulomb potential energy of uniformly charged spheres
@volumes
@njit
def u_coul_ucs(r, rc, z_proj, z_targ):
    q2 = z_proj * z_targ * e2 # MeV.fm
    r_in, r_out = r[r<rc], r[r>=rc] # fm
    uc_in = (3 - r_in*r_in/rc/rc)/rc/2
    uc_out = 1/r_out
    return q2 * append(uc_in, uc_out) # MeV

def u_coul_bifold_d(rho_p_ch_, rho_t_ch_, v_coul_, r, q, R=None, s=None):
    # Direct part of Coulomb potential using folding integrals

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    r_cou = mesh(r[0], 1.25 * r[-1], r[1] - r[0])
    s_cou = mesh(s[0], 1.25 * s[-1], s[1] - s[0])

    s_len = len(s)
    s34 = int(.75 * s_len)
    s_tail = s[s34:]

    v_coul = v_coul_.copy()
    rho_p_ch = rho_p_ch_.copy()
    rho_t_ch = rho_t_ch_.copy()

    # correction for the tail of Coulomb potential
    z_p = rho_p_ch.info[0]['vol2']
    z_t = rho_t_ch.info[0]['vol2']
    u_tail = z_p * z_t * e2 / s_tail


    v_coul.value = interp(s_cou, s, v_coul.value)
    rho_p_ch.value = interp(r_cou, r, rho_p_ch.value)
    rho_t_ch.value = interp(r_cou, r, rho_t_ch.value)

    # actual integration happens here
    u_coul_df_dict = u_bifold_d(rho_p_ch, rho_t_ch, v_coul, r_cou, q, R, s_cou)

    # this is to fix integration errors on the tail part of the coulomb
    # by rescaling folding coulomb to z_p * z_t * e2 * 1/r (analytical coulomb)
    re_scale = u_tail[0]/u_coul_df_dict['func_r']['u_R'][s34]
    u_coul_df_dict['func_r']['u_R'] = append(u_coul_df_dict['func_r']['u_R'][:s34]*re_scale, u_tail)
    u_coul_df_dict['func_i']['u_R'][0]['name'] = 'u_coul_bifold_d'

    return u_coul_df_dict

def vol_msr(R, u_R):
    R2 = R*R
    R4 = R2 * R2
    u_R_vol2 = pi4 * simpson( u_R * R2, R)
    u_R_vol4 = pi4 * simpson( u_R * R4, R)
    u_R_msr  = u_R_vol4 / u_R_vol2
    return u_R_vol2, u_R_vol4, u_R_msr

#..........................................................................#
#******** Density Dependent M3Y - Reid/Paris [B/C/D-DM3Y: Xdm3yn] *********#
#..........................................................................#
def u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul, r, q, R=None, s=None, Cs=1 / 36, dd_name='bdm3y1', vnn_name='reid', u_ex_iter=8):

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    u_d  = u_xdm3yn_d(e_lab, a_proj, rho_p, rho_t, r, q, R, s, dd_name=dd_name, vnn_name=vnn_name)
    u_ex = u_xdm3yn_ex_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_d['func_r']['u_R'], u_coul, r, q, R, s, Cs=Cs,
                          dd_name=dd_name, vnn_name=vnn_name, u_ex_iter=u_ex_iter)
    u_R = u_d['func_r']['u_R'] + u_ex['func_r']['u_R']
    u_R_vol2 = u_d['func_i']['u_R'][0]['vol2'] + u_ex['func_i']['u_R'][0]['vol2']
    u_R_vol4 = u_d['func_i']['u_R'][0]['vol4'] + u_ex['func_i']['u_R'][0]['vol4']
    u_R_msr = u_R_vol4/u_R_vol2
    u_R_info = {'name': f'u_{dd_name}_{vnn_name}_fr', 'L': 0, 'norm': None, 'renorm':1.0,
                'vol2': u_R_vol2, 'vol4': u_R_vol4, 'msr': u_R_msr}

    u_q = pi4 * fourier_with_simpson(u_R, R, q)
    return {'func_i': {'total': {'u_R': [u_R_info]}, 'direct': u_d['func_i'], 'exchange': u_ex['func_i']},
            'func_r': {'total': {'u_R': u_R},        'direct': u_d['func_r'], 'exchange': u_ex['func_r']},
            'func_q': {'total': {'u_R': u_q},        'direct': u_d['func_q'], 'exchange': u_ex['func_q']}}

def u_xdm3yn_zr(e_lab, a_proj, rho_p, rho_t, r, q, R=None, s=None, dd_name='bdm3y1', vnn_name='reid'):

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    u_d  = u_xdm3yn_d(e_lab, a_proj, rho_p, rho_t, r, q, R, s, dd_name=dd_name, vnn_name=vnn_name)
    u_ex = u_xdm3yn_ex_zr(e_lab, a_proj, rho_p, rho_t, r, q, R, s, dd_name=dd_name, vnn_name=vnn_name)
    u_R = u_d['func_r']['u_R'] + u_ex['func_r']['u_R']
    u_R_vol2 = u_d['func_i']['u_R'][0]['vol2'] + u_ex['func_i']['u_R'][0]['vol2']
    u_R_vol4 = u_d['func_i']['u_R'][0]['vol4'] + u_ex['func_i']['u_R'][0]['vol4']
    u_R_msr = u_R_vol4/u_R_vol2
    u_R_info = {'name': f'u_{dd_name}_{vnn_name}_zr', 'L': 0, 'norm': None, 'renorm':1.0,
                'vol2': u_R_vol2, 'vol4': u_R_vol4, 'msr': u_R_msr}

    u_q = pi4 * fourier_with_simpson(u_R, R, q)
    return {'func_i': {'total': {'u_R': [u_R_info]}, 'direct': u_d['func_i'], 'exchange': u_ex['func_i']},
            'func_r': {'total': {'u_R': u_R},        'direct': u_d['func_r'], 'exchange': u_ex['func_r']},
            'func_q': {'total': {'u_R': u_q},        'direct': u_d['func_q'], 'exchange':u_ex['func_q']}}


def u_xdm3yn_d(e_lab, a_proj, rho_p, rho_t, r, q, R=None, s=None, dd_name='bdm3y1', vnn_name='reid'):

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    try:
        c, a, b, g, n, K, implemented = v_xdm3yn_cabgn(dd_name=dd_name, vnn_name=vnn_name)
        if not implemented:
            print(f'{dd_name}_{vnn_name}: (C={c}, alpha={a}, beta={b}, gamma={g}, n={n}) is not implemented yet.')
            quit()
    except KeyError as err:
        print(f'{err} is undefined for {dd_name} - {vnn_name}')
        quit()

    if vnn_name == 'reid':
        gE = 1 - 0.002 * e_lab/a_proj
        vnn = v_m3y_reid_d(s)
    else: # vnn_name == 'paris'
        gE = 1 - 0.003 * e_lab / a_proj
        vnn = v_m3y_paris_d(s)


    u_d_part1 = u_bifold_d(rho_p, rho_t, vnn, r, q, R, s)

    if 'dim3y' in dd_name:
        gE = None
        u_d_part2 = 0
        u_d_part3 = 0
        u_d_part4 = 0
        u_R = u_d_part1['func_r']['u_R']
        u_q = u_d_part1['func_q']['u_R']
    elif 'ddm3y' in dd_name:
        frho_p_dd = f_rho_dd(r, rho_p, beta=b)
        frho_t_dd = f_rho_dd(r, rho_t, beta=b)

        u_d_part2 = u_bifold_d(frho_p_dd, frho_t_dd, vnn, r, q, R, s)
        u_d_part3 = 0
        u_d_part4 = 0
        u_R = c*gE * (u_d_part1['func_r']['u_R'] + a * u_d_part2['func_r']['u_R'])
        u_q = c*gE * (u_d_part1['func_q']['u_R'] + a * u_d_part2['func_q']['u_R'])
    elif 'cdm3y' in dd_name:
        frho_p_dd = f_rho_dd(r, rho_p, beta=b)
        frho_t_dd = f_rho_dd(r, rho_t, beta=b)
        frho_p_bd = f_rho_bd(r, rho_p, n=1)
        frho_t_bd = f_rho_bd(r, rho_t, n=1)

        u_d_part2 = u_bifold_d(frho_p_dd, frho_t_dd, vnn, r, q, R, s)
        u_d_part3 = u_bifold_d(frho_p_bd, rho_t, vnn, r, q, R, s)
        u_d_part4 = u_bifold_d(rho_p, frho_t_bd, vnn, r, q, R, s)
        u_R = c*gE * (u_d_part1['func_r']['u_R'] + a * u_d_part2['func_r']['u_R'] -
                        g*(u_d_part3['func_r']['u_R'] + u_d_part4['func_r']['u_R']))
        u_q = c*gE * (u_d_part1['func_q']['u_R'] + a * u_d_part2['func_q']['u_R'] -
                        g*(u_d_part3['func_q']['u_R'] + u_d_part4['func_q']['u_R']))
    elif 'bdm3y' in dd_name:
        frho_p = f_rho_bd(r, rho_p, n=n)
        frho_t = f_rho_bd(r, rho_t, n=n)

        u_d_part2 = 0
        u_d_part3 = u_bifold_d(frho_p, rho_t, vnn, r, q, R, s)
        u_d_part4 = u_bifold_d(rho_p, frho_t, vnn, r, q, R, s)
        u_R = c*gE * (u_d_part1['func_r']['u_R'] - g * (u_d_part3['func_r']['u_R'] + u_d_part4['func_r']['u_R']))
        u_q = c*gE * (u_d_part1['func_q']['u_R'] - g * (u_d_part3['func_q']['u_R'] + u_d_part4['func_q']['u_R']))
    else:
        print('Something went wrong!')
        quit()

    u_R_vol2, u_R_vol4, u_R_msr = vol_msr(R, u_R)
    u_R_info = {'name': f'u_{dd_name}_{vnn_name}_d', 'L': 0, 'norm': None, 'renorm':1.0,
                'vol2': u_R_vol2, 'vol4': u_R_vol4, 'msr': u_R_msr,
                'c':c, 'alpha':a, 'beta':b, 'gamma':g, 'n':n, 'gE':gE}

    if 'dim3y' in dd_name:
        return {'func_i': {'u_R': [u_R_info], 'part1': u_d_part1['func_i']},
                'func_r': {'u_R': u_R,        'part1': u_d_part1['func_r']},
                'func_q': {'u_R': u_q,        'part1': u_d_part1['func_q']}}
    elif 'ddm3y' in dd_name:
        return {'func_i': {'u_R': [u_R_info], 'part1': u_d_part1['func_i'], 'part2': u_d_part2['func_i']},
                'func_r': {'u_R': u_R,        'part1': u_d_part1['func_r'], 'part2': u_d_part2['func_r']},
                'func_q': {'u_R': u_q,        'part1': u_d_part1['func_q'], 'part2': u_d_part2['func_q']}}
    elif 'cdm3y' in dd_name:
        return {'func_i': {'u_R': [u_R_info], 'part1': u_d_part1['func_i'], 'part2': u_d_part2['func_i'], 'part3': u_d_part3['func_i'], 'part4': u_d_part4['func_i']},
                'func_r': {'u_R': u_R,        'part1': u_d_part1['func_r'], 'part2': u_d_part2['func_r'], 'part3': u_d_part3['func_r'], 'part4': u_d_part4['func_r']},
                'func_q': {'u_R': u_q,        'part1': u_d_part1['func_q'], 'part2': u_d_part2['func_q'], 'part3': u_d_part3['func_q'], 'part4': u_d_part4['func_q']}}
    elif 'bdm3y' in dd_name:
        return {'func_i': {'u_R': [u_R_info], 'part1': u_d_part1['func_i'], 'part3': u_d_part3['func_i'], 'part4': u_d_part4['func_i']},
                'func_r': {'u_R': u_R,        'part1': u_d_part1['func_r'], 'part3': u_d_part3['func_r'], 'part4': u_d_part4['func_r']},
                'func_q': {'u_R': u_q,        'part1': u_d_part1['func_q'], 'part3': u_d_part3['func_q'], 'part4': u_d_part4['func_q']}}


def u_xdm3yn_ex_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_d, u_coul_dict, r, q, R=None, s=None,
                   Cs=1/36, dd_name='bdm3y1', vnn_name='reid', u_ex_iter=8):

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    # folded charge densities to obtain Coulomb
    # or uniformly charged spheres to obtain Coulomb
    try:
        u_coul = u_coul_dict()
        u_coul_info = u_coul_dict.info
    except TypeError:
        u_coul = u_coul_dict['func_r']['u_R']
        u_coul_info = u_coul_dict['func_i']#['u_R']
    except:
        print("Please use 'u_coul_ucs()' or 'u_coul_bifold_d()' for Coulomb potential.")

    try:
        c, a, b, g, n, K, implemented = v_xdm3yn_cabgn(dd_name=dd_name, vnn_name=vnn_name)
        if not implemented:
            print(f'{dd_name}_{vnn_name}: (C={c}, alpha={a}, beta={b}, gamma={g}, n={n}) is not implemented yet.')
            quit()
    except KeyError as err:
        print(f'{err} is undefined for {dd_name} - {vnn_name}')
        quit()

    if vnn_name == 'reid':
        gE = 1 - 0.002 * e_lab/a_proj
        vnn = v_m3y_reid_ex_fr(s)
    else: # vnn_name == 'paris'
        gE = 1 - 0.003 * e_lab / a_proj
        vnn = v_m3y_paris_ex_fr(s)

    a_total = a_proj + a_targ
    ecm = e_lab * a_targ / a_total
    a_reduced = a_targ * a_proj / a_total

    rho_p_name = rho_p.info[0]['name']
    rho_t_name = rho_t.info[0]['name']
    check_f_names = ['f_external', 'f_internet', 'f_nudat']

    kf_p = k_fermi_spline(r, rho_p(), Cs=Cs) if any([ci in rho_p_name for ci in check_f_names]) else k_fermi(r, rho_p(), Cs=Cs)
    kf_t = k_fermi_spline(r, rho_t(), Cs=Cs) if any([ci in rho_t_name for ci in check_f_names]) else k_fermi(r, rho_t(), Cs=Cs)

    fa = pi4 * fqs_with_simpson(rho_p(), r, kf_p(), s, q)
    fA = pi4 * fqs_with_simpson(rho_t(), r, kf_t(), s, q)

    if 'dim3y' in dd_name:
        gE = None
        dFqs = fa * fA
    elif 'ddm3y' in dd_name:
        frho_p_dd = f_rho_dd(r, rho_p, beta=b)
        frho_t_dd = f_rho_dd(r, rho_t, beta=b)

        fa_exp = pi4 * fqs_with_simpson(frho_p_dd(), r, kf_p(), s, q)
        fA_exp = pi4 * fqs_with_simpson(frho_t_dd(), r, kf_t(), s, q)

        dFqs = fa * fA + a * (fa_exp * fA_exp)
    elif 'cdm3y' in dd_name:
        frho_p_dd = f_rho_dd(r, rho_p, beta=b)
        frho_t_dd = f_rho_dd(r, rho_t, beta=b)
        frho_p_bd = f_rho_bd(r, rho_p, n=1)
        frho_t_bd = f_rho_bd(r, rho_t, n=1)
        fa_exp = pi4 * fqs_with_simpson(frho_p_dd(), r, kf_p(), s, q)
        fA_exp = pi4 * fqs_with_simpson(frho_t_dd(), r, kf_t(), s, q)
        fa2 = pi4 * fqs_with_simpson(frho_p_bd(), r, kf_p(), s, q)
        fA2 = pi4 * fqs_with_simpson(frho_t_bd(), r, kf_t(), s, q)
        dFqs = fa * fA + a * (fa_exp * fA_exp) - g* (fa * fA2 + fa2 * fA)
    elif 'bdm3y' in dd_name:
        frho_p_bd = f_rho_bd(r, rho_p, n=n)
        frho_t_bd = f_rho_bd(r, rho_t, n=n)
        fa2 = pi4 * fqs_with_simpson(frho_p_bd(), r, kf_p(), s, q)
        fA2 = pi4 * fqs_with_simpson(frho_t_bd(), r, kf_t(), s, q)

        dFqs = fa * fA - g * (fa * fA2 + fa2 * fA)
    else:
        print('Something went wrong!')
        quit()

    dGRs = pi2_inv * gRs_with_simpson(dFqs, R, s, q)


    if u_ex_iter <= 0:
        print(f'u_ex_iter must be at least 1! It is raised from {u_ex_iter} to 1.')
        u_ex_iter = 1

    u_nuc = u_coul + u_d
    u_ex = 0
    for i in range(u_ex_iter):
        k2_local_mom_direct = 2 * mu_c2 * a_reduced / hbc / hbc * (ecm - (u_nuc + u_ex))
        # abs is not exist in original definition
        # it is here for keep k_local real valued!
        k_local = sqrt( abs(k2_local_mom_direct) ) / a_reduced
        if c==None:
            # density independent finite range exchange potential
            u_ex = pi4 * u_ex_with_simpson(dGRs, k_local, vnn(), R, s)
        else:
            u_ex = c * gE * pi4 * u_ex_with_simpson(dGRs, k_local, vnn(), R, s)

    u_R = u_ex
    u_R_vol2, u_R_vol4, u_R_msr = vol_msr(R, u_R)
    u_R_info = {'name': f'u_{dd_name}_{vnn_name}_ex_fr', 'L': 0, 'norm': None, 'renorm':1.0,
                'vol2': u_R_vol2, 'vol4': u_R_vol4, 'msr': u_R_msr,
                'c':c, 'alpha':a, 'beta':b, 'gamma':g, 'n':n, 'gE':gE}

    u_q = pi4 * fourier_with_simpson(u_R, R, q)
    rho_pq = pi4 * fourier_with_simpson(rho_p(), r, q)
    rho_tq = pi4 * fourier_with_simpson(rho_t(), r, q)
    vnn_q  = pi4 * fourier_with_simpson(vnn(), s, q)
    kf_pq = pi4 * fourier_with_simpson(kf_p(), r, q)
    kf_tq = pi4 * fourier_with_simpson(kf_t(), r, q)

    return {'func_i': {'u_R': [u_R_info], 'rho_p':rho_p.info, 'rho_t':rho_t.info,
                       'vnn':vnn.info,    'kf_p': kf_p.info,  'kf_t': kf_t.info, 'u_coul': u_coul_info},
            'func_r': {'u_R': u_R,        'rho_p':rho_p(),    'rho_t':rho_t(),
                       'vnn':vnn(),       'kf_p': kf_p(),     'kf_t': kf_t(),    'u_coul': u_coul},
            'func_q': {'u_R': u_q, 'rho_p': rho_pq, 'rho_t': rho_tq, 'vnn': vnn_q, 'kf_p': kf_pq, 'kf_t': kf_tq}}

def u_xdm3yn_ex_zr(e_lab, a_proj, rho_p, rho_t, r, q, R=None, s=None, dd_name='bdm3y1', vnn_name='reid'):

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    try:
        c, a, b, g, n, K, implemented = v_xdm3yn_cabgn(dd_name=dd_name, vnn_name=vnn_name)
        if not implemented:
            print(f'{dd_name}_{vnn_name}: (C={c}, alpha={a}, beta={b}, gamma={g}, n={n}) is not implemented yet.')
            quit()
    except KeyError as err:
        print(f'{err} is undefined for {dd_name} - {vnn_name}')
        quit()

    if vnn_name == 'reid':
        gE = 1 - 0.002 * e_lab/a_proj
        vnn = v_m3y_reid_ex_zr(s, e_lab, a_proj, L=0)
    else: # vnn_name == 'paris'
        gE = 1 - 0.003 * e_lab / a_proj
        vnn = v_m3y_paris_ex_zr(s, e_lab, a_proj, L=0)


    u_d_part1 = u_bifold_ex_zr(rho_p, rho_t, vnn, r, q, R, s)

    if 'dim3y' in dd_name:
        gE = None
        u_d_part2 = 0
        u_d_part3 = 0
        u_d_part4 = 0
        u_R = u_d_part1['func_r']['u_R']
        u_q = u_d_part1['func_q']['u_R']
    elif 'ddm3y' in dd_name:
        frho_p_dd = f_rho_dd(r, rho_p, beta=b)
        frho_t_dd = f_rho_dd(r, rho_t, beta=b)
        u_d_part2 = u_bifold_ex_zr(frho_p_dd, frho_t_dd, vnn, r, q, R, s)
        u_d_part3 = 0
        u_d_part4 = 0
        u_R = c*gE * (u_d_part1['func_r']['u_R'] + a * u_d_part2['func_r']['u_R'])
        u_q = c*gE * (u_d_part1['func_q']['u_R'] + a * u_d_part2['func_q']['u_R'])
    elif 'cdm3y' in dd_name:
        frho_p_dd = f_rho_dd(r, rho_p, beta=b)
        frho_t_dd = f_rho_dd(r, rho_t, beta=b)
        frho_p_bd = f_rho_bd(r, rho_p, n=1)
        frho_t_bd = f_rho_bd(r, rho_t, n=1)

        u_d_part2 = u_bifold_ex_zr(frho_p_dd, frho_t_dd, vnn, r, q, R, s)
        u_d_part3 = u_bifold_ex_zr(frho_p_bd, rho_t, vnn, r, q, R, s)
        u_d_part4 = u_bifold_ex_zr(rho_p, frho_t_bd, vnn, r, q, R, s)
        u_R = c*gE * (u_d_part1['func_r']['u_R'] + a * u_d_part2['func_r']['u_R'] -
                        g*(u_d_part3['func_r']['u_R'] + u_d_part4['func_r']['u_R']))
        u_q = c*gE * (u_d_part1['func_q']['u_R'] + a * u_d_part2['func_q']['u_R'] -
                        g*(u_d_part3['func_q']['u_R'] + u_d_part4['func_q']['u_R']))
    elif 'bdm3y' in dd_name:
        frho_p = f_rho_bd(r, rho_p, n=n)
        frho_t = f_rho_bd(r, rho_t, n=n)
        u_d_part2 = 0
        u_d_part3 = u_bifold_ex_zr(frho_p, rho_t, vnn, r, q, R, s)
        u_d_part4 = u_bifold_ex_zr(rho_p, frho_t, vnn, r, q, R, s)
        u_R = c*gE * (u_d_part1['func_r']['u_R'] - g * (u_d_part3['func_r']['u_R'] + u_d_part4['func_r']['u_R']))
        u_q = c*gE * (u_d_part1['func_q']['u_R'] - g * (u_d_part3['func_q']['u_R'] + u_d_part4['func_q']['u_R']))
    else:
        print('Something went wrong!')
        quit()

    u_R_vol2, u_R_vol4, u_R_msr = vol_msr(R, u_R)
    u_R_info = {'name': f'u_{dd_name}_{vnn_name}_ex_zr', 'L': 0, 'norm': None, 'renorm':1.0,
                'vol2': u_R_vol2, 'vol4': u_R_vol4, 'msr': u_R_msr,
                'c':c, 'alpha':a, 'beta':b, 'gamma':g, 'n':n, 'gE':gE}
    if 'dim3y' in dd_name:
        return {'func_i': {'u_R': [u_R_info], 'part1': u_d_part1['func_i']},
                'func_r': {'u_R': u_R,        'part1': u_d_part1['func_r']},
                'func_q': {'u_R': u_q,        'part1': u_d_part1['func_q']}}
    elif 'ddm3y' in dd_name:
        return {'func_i': {'u_R': [u_R_info], 'part1': u_d_part1['func_i'], 'part2': u_d_part2['func_i']},
                'func_r': {'u_R': u_R,        'part1': u_d_part1['func_r'], 'part2': u_d_part2['func_r']},
                'func_q': {'u_R': u_q,        'part1': u_d_part1['func_q'], 'part2': u_d_part2['func_q']}}
    elif 'cdm3y' in dd_name:
        return {'func_i': {'u_R': [u_R_info], 'part1': u_d_part1['func_i'], 'part2': u_d_part2['func_i'], 'part3': u_d_part3['func_i'], 'part4': u_d_part4['func_i']},
                'func_r': {'u_R': u_R,        'part1': u_d_part1['func_r'], 'part2': u_d_part2['func_r'], 'part3': u_d_part3['func_r'], 'part4': u_d_part4['func_r']},
                'func_q': {'u_R': u_q,        'part1': u_d_part1['func_q'], 'part2': u_d_part2['func_q'], 'part3': u_d_part3['func_q'], 'part4': u_d_part4['func_q']}}
    elif 'bdm3y' in dd_name:
        return {'func_i': {'u_R': [u_R_info], 'part1': u_d_part1['func_i'], 'part3': u_d_part3['func_i'], 'part4': u_d_part4['func_i']},
                'func_r': {'u_R': u_R,        'part1': u_d_part1['func_r'], 'part3': u_d_part3['func_r'], 'part4': u_d_part4['func_r']},
                'func_q': {'u_R': u_q,        'part1': u_d_part1['func_q'], 'part3': u_d_part3['func_q'], 'part4': u_d_part4['func_q']}}


#### LOCAL FERMI MOMENTUM
@volumes
@njit
def k_fermi(r, f, Cs=1/36):
    """Calculates local momentum of 'f' function.
    """
    df = f_der1(r, f)
    ddf = f_der2(r, f)

    kf1 = power(3*pi_sqr/2 * f, 2/3)
    kf2_num = 5*Cs*power(df, 2)
    kf2_denom = 3*power(f, 2)
    kf3_num = 5*ddf
    kf3_denom = 36*f

    kf_sqr = kf1 + kf2_num/kf2_denom +  kf3_num/kf3_denom
    kf = sqrt(kf_sqr)

    kf[-1] = f__1(r, kf)

    return kf

@volumes
def k_fermi_spline(r, f, Cs=1/36):
    """Calculates local momentum when 'f' functions is one oft the 'f_external', 'f_internet' and 'f_nudat' functions.
    Spline is necessary to smooth the derivatives of the 'f' function.
    """
    f = f_spline(r, f)
    df = f_der1s(r, f)
    ddf = f_der2s(r, f)

    kf1 = f_spline(r, power(3*pi_sqr/2 * f, 2/3))
    kf2_num = f_spline(r, 5*Cs*power(df, 2))
    kf2_denom = f_spline(r, 3*power(f, 2))
    kf3_num = 5*ddf
    kf3_denom = 36*f

    kf_sqr = f_spline(r, kf1 + kf2_num/kf2_denom + kf3_num/kf3_denom)
    kf = f_spline(r, sqrt(kf_sqr))

    return kf
#..........................................................................#
#******** Density Dependent M3Y - Reid/Paris [B/C/D-DM3Y: Xdm3yn] *********#
#..........................................................................#


#..........................................................................#
#****************** Density Dependent M3Y - Reid [DDM3Y] ******************#
#..........................................................................#
def u_ddm3y_reid_zr(e_lab, a_proj, rho_p, rho_t, r, q, R=None, s=None):

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    u_d = u_ddm3y_reid_d(e_lab, a_proj, rho_p, rho_t, r, q, R, s)
    u_ex = u_ddm3y_reid_ex_zr(e_lab, a_proj, rho_p, rho_t, r, q, R, s)
    u_R = u_d['func_r']['u_R'] + u_ex['func_r']['u_R']
    u_q = u_d['func_q']['u_R'] + u_ex['func_q']['u_R']
    u_R_vol2 = u_d['func_i']['u_R'][0]['vol2'] + u_ex['func_i']['u_R'][0]['vol2']
    u_R_vol4 = u_d['func_i']['u_R'][0]['vol4'] + u_ex['func_i']['u_R'][0]['vol4']
    u_R_msr = u_R_vol4/u_R_vol2
    u_R_info = {'name': 'u_ddm3y_reid_zr', 'L': 0, 'norm': None, 'renorm': 1.0,
                'vol2': u_R_vol2, 'vol4': u_R_vol4, 'msr': u_R_msr}

    return {'func_i': {'total': {'u_R': [u_R_info]}, 'direct': u_d['func_i'], 'exchange': u_ex['func_i']},
            'func_r': {'total': {'u_R': u_R},        'direct': u_d['func_r'], 'exchange': u_ex['func_r']},
            'func_q': {'total': {'u_R': u_q},        'direct': u_d['func_q'], 'exchange': u_ex['func_q']}}


def u_ddm3y_reid_d(e_lab, a_proj, rho_p, rho_t, r, q, R=None, s=None):

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    # DDM3Y zero range direct part
    # DDM3Y energy dependent parameters
    c, a, b = v_ddm3y_reid_cab(e_lab, a_proj)

    frho_p = f_rho_dd(r, rho_p, beta=b)
    frho_t = f_rho_dd(r, rho_t, beta=b)
    vnn = v_m3y_reid_d(r)

    u_d_part1 = u_bifold_d(rho_p, rho_t, vnn, r, q, R, s)
    u_d_part2 = u_bifold_d(frho_p, frho_t, vnn, r, q, R, s)

    u_R = c*u_d_part1['func_r']['u_R'] + c*a*u_d_part2['func_r']['u_R']
    u_q = c*u_d_part1['func_q']['u_R'] + c*a*u_d_part2['func_q']['u_R']
    u_R_vol2, u_R_vol4, u_R_msr = vol_msr(R, u_R)
    u_R_info = {'name': 'u_ddm3y_reid_d', 'L': 0, 'norm': None, 'renorm':1.0,
                'vol2': u_R_vol2, 'vol4': u_R_vol4, 'msr': u_R_msr,
                'c':c, 'alpha':a, 'beta':b}

    return {'func_i': {'u_R': [u_R_info], 'part1': u_d_part1['func_i'], 'part2': u_d_part2['func_i']},
            'func_r': {'u_R': u_R,        'part1': u_d_part1['func_r'], 'part2': u_d_part2['func_r']},
            'func_q': {'u_R': u_q,        'part1': u_d_part1['func_q'], 'part2': u_d_part2['func_q']}}


def u_ddm3y_reid_ex_zr(e_lab, a_proj, rho_p, rho_t, r, q, R=None, s=None):

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    # DDM3Y zero range exchange part
    # DDM3Y energy dependent parameters
    c, a, b = v_ddm3y_reid_cab(e_lab, a_proj)

    frho_p = f_rho_dd(r, rho_p, beta=b)
    frho_t = f_rho_dd(r, rho_t, beta=b)
    vnn = v_m3y_reid_ex_zr(r,e_lab, a_proj)

    u_e_part1 = u_bifold_ex_zr(rho_p, rho_t, vnn, r, q, R, s)
    u_e_part2 = u_bifold_ex_zr(frho_p, frho_t, vnn, r, q, R, s)

    u_R = c*u_e_part1['func_r']['u_R'] + c*a*u_e_part2['func_r']['u_R']
    u_q = c*u_e_part1['func_q']['u_R'] + c*a*u_e_part2['func_q']['u_R']
    u_R_vol2, u_R_vol4, u_R_msr = vol_msr(R, u_R)
    u_R_info = {'name': 'u_ddm3y_reid_ex_zr', 'L': 0, 'norm': None, 'renorm':1.0,
                'vol2': u_R_vol2, 'vol4': u_R_vol4, 'msr': u_R_msr,
                'c':c, 'alpha':a, 'beta':b}
    return {'func_i': {'u_R': [u_R_info], 'part1': u_e_part1['func_i'], 'part2': u_e_part2['func_i']},
            'func_r': {'u_R': u_R,        'part1': u_e_part1['func_r'], 'part2': u_e_part2['func_r']},
            'func_q': {'u_R': u_q,        'part1': u_e_part1['func_q'], 'part2': u_e_part2['func_q']}}
#..........................................................................#
#****************** Density Dependent M3Y - Reid [DDM3Y] ******************#
#..........................................................................#


#..........................................................................#
#****************** Density Independent M3Y - Reid/Paris ******************#
#..........................................................................#
def u_m3y_reid_zr(e_lab, a_proj, rho_p, rho_t, r, q, R=None, s=None):
    
    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    vnn_d = v_m3y_reid_d(r)
    vnn_ex = v_m3y_reid_ex_zr(r, e_lab, a_proj, L=0)
    u_m3y_zr_dict = u_bifold_zr(rho_p, rho_t, vnn_d, vnn_ex, r, q, R, s)
    u_m3y_zr_dict['func_i']['total']['u_R'][0]['name'] = 'u_m3y_reid_zr'
    return u_m3y_zr_dict

def u_m3y_paris_zr(e_lab, a_proj, rho_p, rho_t, r, q, R=None, s=None):

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    vnn_d = v_m3y_paris_d(r)
    vnn_ex = v_m3y_paris_ex_zr(r, e_lab, a_proj, L=0)
    u_m3y_zr_dict = u_bifold_zr(rho_p, rho_t, vnn_d, vnn_ex, r, q, R, s)
    u_m3y_zr_dict['func_i']['total']['u_R'][0]['name'] = 'u_m3y_paris_zr'
    return u_m3y_zr_dict


def u_bifold_zr(rho_p, rho_t, vnn_d, vnn_ex, r, q, R=None, s=None):

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    u_d_dict = u_bifold_d(rho_p, rho_t, vnn_d, r, q, R, s)
    u_ex_zr_dict = u_bifold_ex_zr(rho_p, rho_t, vnn_ex, r, q, R, s)
    u_R = u_d_dict['func_r']['u_R'] + u_ex_zr_dict['func_r']['u_R']
    u_q = u_d_dict['func_q']['u_R'] + u_ex_zr_dict['func_q']['u_R']
    u_R_vol2 = u_d_dict['func_i']['u_R'][0]['vol2'] + u_ex_zr_dict['func_i']['u_R'][0]['vol2']
    u_R_vol4 = u_d_dict['func_i']['u_R'][0]['vol4'] + u_ex_zr_dict['func_i']['u_R'][0]['vol4']
    u_R_msr = u_R_vol4/u_R_vol2
    u_R_info = {'name':'u_bifold_zr', 'L':0, 'norm': None, 'renorm':1.0,
                'vol2':u_R_vol2, 'vol4':u_R_vol4, 'msr':u_R_msr}
    return {'func_i': {'total': {'u_R': [u_R_info]}, 'direct': u_d_dict['func_i'], 'exchange': u_ex_zr_dict['func_i']},
            'func_r': {'total': {'u_R': u_R},        'direct': u_d_dict['func_r'], 'exchange': u_ex_zr_dict['func_r']},
            'func_q': {'total': {'u_R': u_q},        'direct': u_d_dict['func_q'], 'exchange': u_ex_zr_dict['func_q']}}

def u_bifold_d(rho_p, rho_t, vnn, r, q, R=None, s=None):

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    rho_pq = pi4 * fourier_with_simpson(rho_p(), r, q)
    rho_tq = pi4 * fourier_with_simpson(rho_t(), r, q)
    vnn_q  = pi4 * fourier_with_simpson(vnn(), s, q)

    u_q = rho_pq * rho_tq * vnn_q
    u_R = pi2_inv * fourier_with_simpson(u_q, q, R)
    u_R_vol2, u_R_vol4, u_R_msr = vol_msr(R, u_R)

    u_R_info = {'name':'u_direct', 'L':0, 'norm': None, 'renorm':1.0,
                'vol2':u_R_vol2, 'vol4':u_R_vol4, 'msr':u_R_msr}
    return {'func_i': {'u_R': [u_R_info], 'rho_p': rho_p.info, 'rho_t': rho_t.info, 'vnn': vnn.info},
            'func_r': {'u_R': u_R,        'rho_p': rho_p(),    'rho_t': rho_t(),    'vnn': vnn()},
            'func_q': {'u_R': u_q,        'rho_p': rho_pq,     'rho_t': rho_tq,     'vnn': vnn_q}}


def u_bifold_ex_zr(rho_p, rho_t, vnn, r, q, R=None, s=None):

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    rho_pq = pi4 * fourier_with_simpson(rho_p(), r, q)
    rho_tq = pi4 * fourier_with_simpson(rho_t(), r, q)
    vnn_q  = vnn()[0] + 0*q

    u_q = rho_pq * rho_tq * vnn_q
    u_R = pi2_inv * fourier_with_simpson(u_q, q, R)
    u_R_vol2, u_R_vol4, u_R_msr = vol_msr(R, u_R)

    u_R_info = {'name':'u_exchange_zr', 'L':0, 'norm': None, 'renorm':1.0,
                'vol2':u_R_vol2, 'vol4':u_R_vol4, 'msr':u_R_msr}
    return {'func_i': {'u_R': [u_R_info], 'rho_p': rho_p.info, 'rho_t': rho_t.info, 'vnn': vnn.info},
            'func_r': {'u_R': u_R,        'rho_p': rho_p(),    'rho_t': rho_t(),    'vnn': vnn()},
            'func_q': {'u_R': u_q,        'rho_p': rho_pq,     'rho_t': rho_tq,     'vnn': vnn_q}}
#..........................................................................#
#****************** Density Independent M3Y - Reid/Paris ******************#
#..........................................................................#
