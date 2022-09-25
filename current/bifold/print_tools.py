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
This module contains the printing tools for BiFold.
"""


DASHED = '- ' * 47
SOLID = '─' * 50
DOUBLE = '═' * 50
TABLE_HEADER = 'density/interaction               L          norm        renorm      vol2        vol4            msr'
TABLE_UNDERLINE =  '─' * 110

def print_title(title='', n_l='\n', omit=''):
    l_title = len(title)
    sides = (110-l_title)//2
    if sides>1:
        side_strip = '═' * (sides-1)
        print(f'{omit}{n_l}{omit}{side_strip}╣{title}╠{side_strip}{omit}{n_l}{omit}')
    else:
        print(f'{omit}{n_l}{omit}{title}')

def print_section(val, l_right=SOLID, l_fet=SOLID, n_l='\n', omit=''):
    print(f'{omit}{n_l}{omit}{l_fet} {val} {l_right}')

def print_subsection(val, fmts='10s', l_left='───────', l_right=DASHED, n_l='\n', omit=''):
    print(f'{omit}{n_l}{omit}{l_left}{val:{fmts}}{l_right}')


def print_bifold_logo(omit=''):
    # https://textkool.com/en/ascii-art-generator?hl=default&vl=default&font=Doh&text=BF
    name_bifold = \
"""
BBBBBBBBBBBBBBBBB  FFFFFFFFFFFFFFFFFFFFFF
B                B F                    F
B      BBBBBB     BF                    F
BB     B     B     FF      FFFFFFFFF    F
  B    B     B     B F     F       FFFFFF
  B    B     B     B F     F             
  B    BBBBBB     B  F      FFFFFFFFFF   
  B             BB   F               F   
  B    BBBBBB     B  F               F   
  B    B     B     B F      FFFFFFFFFF   
  B    B     B     B F     F             
  B    B     B     B F     F             
BB     BBBBBB      FF       FF           
B                 BF        FF           
B                B F        FF           
BBBBBBBBBBBBBBBBB  FFFFFFFFFFF           """

    for fn in name_bifold.split('\n'):
        print(f'{omit}{fn:^110s}')


def print_ncol(r, f, ncol=4, show='u_R_short', fmt=None, header=None, underline = '─', omit='', data_format=0):
    '''
    It prints r vs f as a tabular.

    show == 'u_R_short' show the folded potential (summarized) v.s. r including 'info'
    show == 'all_short' show everything calculated in a summarized format
    '''

    fmt = ['>8.3f', '>17.9f'] if fmt == None else fmt
    header = ['r', 'f(r)'] if header == None else header

    if data_format==0:
        table = list(map(lambda x: f'{x[0]:{fmt[0]}} {x[1]:{fmt[1]}}\t', zip(r, f)))
    else:
        table = list(map(lambda x: f'{x:{fmt[1]}}\t', f))

    hr = len(f'{r[0]:{fmt[0]}}')
    hf = len(f'{f[0]:{fmt[1]}}')

    if data_format==0:
        header_fmt = f'{header[0]:^{hr}s} {header[1]:^{hf}s}\t'
        header_u = f" {underline * (hr - 1)}   {underline * (hf - 2)}\t"*ncol
    else:
        header_fmt = f'{header[1]:^{hf}s}\t'
        header_u = f" {underline * (hf - 2)}\t"*ncol

    header_str = ''.join([header_fmt for i in range(ncol)])

    nrow = len(table) // ncol
    dr = r[1] - r[0]
    n = nrow * ncol
    s_dr = f'd{header[0]}'
    s_r0 = f'{header[0]}0'
    if data_format==1:
        print(f'{omit}{s_dr:^{hf}s}\n{omit}{s_r0:^{hf}s}')
        print(f'{dr:{fmt[0]}}\n{r[0]:{fmt[0]}}')
    elif data_format==2:
        print(f'{omit}{"n":^{hf}s}\n{omit}{s_dr:^{hf}s}\n{omit}{s_r0:^{hf}s}')
        print(f'{n:{hr}d}\n{dr:{fmt[0]}}\n{r[0]:{fmt[0]}}')

    print(f'{omit}\n{omit}{header_str[(len(omit)):]}')
    print(f'{omit}{header_u[len(omit):]}')

    if show[-6:]=='_short':
        for j, z in enumerate(zip(*(table[nrow * i:nrow * (i + 1)] for i in range(ncol)))):
            if j<4 or j>nrow-5:
                print(''.join(z))
            elif 4<j<8:
                print(f"{'...':>{hr}s}")
    else:
        for z in zip(*(table[nrow * i:nrow * (i + 1)] for i in range(ncol))):
            print(''.join(z))


def print_ncol_loop(k, v, r, R, s, show='u_R_short', ncol=4, fmt=None, omit='', data_format=0):
    '''
    It prints r vs f as a tabular using print_ncol.
    This is necessary to be able to print all the components of a folded potential.

    show == 'u_R_short' show the folded potential (summarized) v.s. r including 'info'
    show == 'all_short' show everything calculated in a summarized format
    '''
    if k=='u_R':
        print_ncol(R, v, header=['R', k], show=show, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)
        if show == 'u_R' or show == 'u_R_short':
            return
    elif 'u_coul' in k:
        print_ncol(R, v, header=['R', k], show=show, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)
    elif k == 'rho_p' or k == 'rho_t' or k=='kf_p' or k=='kf_t':
        print_ncol(r, v, header=['r', k], show=show, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)
    elif k == 'vnn':
        print_ncol(s, v, header=['s', k], show=show, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)
    else:
        print(f'Something went wrong while printing! {k} is undefined')
        return()


def print_func_info(key, val):
    '''
    It prints info of a folded potential: names, norm, renorm, vol2, vol4, msr
    '''
    for vali in val:
        line_length = 0
        for ki, vi in vali.items():
            if line_length < len(TABLE_UNDERLINE) - 8:
                if isinstance(vi, str):
                    line = f'{key:9s}: {vi:20s}'
                    print(line, end=' ')
                elif isinstance(vi, int):
                    line = f'{vi:3d}'
                    print(line, end=' ')
                elif isinstance(vi, float):
                    line = f'{vi:12.3f}'
                    print(line, end=' ')
                else:
                    line = f'{str(vi):>13s}'
                    print(line, end=' ')
                line_length += len(line) + 1
            else:
                if isinstance(vi, str):
                    line = f'{ki:>5s} : {vi:<20s}'
                    print(line, end=' ')
                elif isinstance(vi, int):
                    line = f'{ki:>5s} = {vi:<3d}'
                    print(line, end=' ')
                elif isinstance(vi, float):
                    line = f'{ki:>5s} = {vi:<8.3f}'
                    print(line, end=' ')
                elif vi!=None:
                    line = f'{ki:>5s} : {str(vi):<13s}'
                    print(line, end=' ')
        print()


def print_bifold(u, r, q, R=None, s=None, show='u_R_short', info=True, ncol=4, fmt=None, omit='', data_format=0):
    '''
    It prints only direct/exchange part of the density independent folded potentials.

    if info == True then
                    all definitions below is True
                    else
                    all definitions below is True with no 'info' printed.
    show == 'info' only show names, norm, renorm, vol2, vol4, msr
    show == 'u_R' show the folded potential v.s. r including 'info'
    show == 'u_R_short' show the folded potential (summarized) v.s. r including 'info'
    show == 'all' show everything calculated
    show == 'all_short' show everything calculated in a summarized format
    '''

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    if info:
        print(TABLE_HEADER)
        print( TABLE_UNDERLINE)
        for k, v in u['func_i'].items():
            print_func_info(k, v)

    if show == 'info':
        print(omit)
        return

    if not (show == 'u_R' or show == 'u_R_short'):
        print_section('r  space', omit=omit)
    for k, v in u['func_r'].items():
        print_ncol_loop(k, v, r, R, s, show=show, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)
        if show == 'u_R' or show == 'u_R_short':
            print(omit)
            return

    print_section('q  space', omit=omit)
    for k, v in u['func_q'].items():
        print_ncol(q, v, header=['q', k], show=show, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)


def print_bifold_di(u, r, q, R=None, s=None, show='u_R_short', info=True, ncol=4, fmt=None, omit='', data_format=0):
    '''
    It prints direct and exchange parts of the density independent folded potentials.

    show == 'info' only show names, norm, renorm, vol2, vol4, msr
    show == 'u_R' show the folded potential v.s. r including 'info'
    show == 'u_R_short' show the folded potential (summarized) v.s. r including 'info'
    show == 'all' show everything calculated
    show == 'all_short' show everything calculated in a summarized format
    '''

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    if info:
        print(TABLE_HEADER)
        print(TABLE_UNDERLINE)
        for ki, vi in u['func_i'].items():
            print_subsection(ki, n_l='', omit=omit)
            for kj, vj in vi.items():
                print_func_info(kj, vj)

    if show == 'info':
        print(omit)
        return

    if not (show == 'u_R' or show == 'u_R_short'):
        print_section('r  space', omit=omit)
    for ki, vi in u['func_r'].items():
        if not (show == 'u_R' or show == 'u_R_short'):
            print_subsection(ki, omit=omit)
        for kj, vj in vi.items():
            print_ncol_loop(kj, vj, r, R, s, show=show, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)
            if show == 'u_R' or show == 'u_R_short':
                print(omit)
                return

    print_section('q  space', omit=omit)
    for ki, vi in u['func_q'].items():
        print_subsection(ki, omit=omit)
        for kj, vj in vi.items():
            print_ncol(q, vj, header=['q', kj], show=show, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)

    print(omit)


def print_bifold_dd(u, r, q, R=None, s=None, show='u_R_short', info=True, ncol=4, fmt=None, omit='', data_format=0):
    '''
    It prints direct and exchange parts of the density dependent folded potentials.

    show == 'info' only show names, norm, renorm, vol2, vol4, msr
    show == 'u_R' show the folded potential v.s. r including 'info'
    show == 'u_R_short' show the folded potential (summarized) v.s. r including 'info'
    show == 'all' show everything calculated
    show == 'all_short' show everything calculated in a summarized format
    '''

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    if info:
        print(TABLE_HEADER)
        print( TABLE_UNDERLINE)
        for ki, vi in u['func_i'].items():
            print_subsection(ki, n_l='', omit=omit)
            for kj, vj in vi.items():
                try:
                    print_func_info(kj, vj)
                except AttributeError:
                    print_subsection(kj, n_l='', omit=omit)
                    for km, vm in vj.items():
                        print_func_info(km, vm)


    if show == 'info':
        print(omit)
        return

    if not (show == 'u_R' or show == 'u_R_short'):
        print_section('r  space', omit=omit)
    for ki, vi in u['func_r'].items():
        if not (show == 'u_R' or show == 'u_R_short'):
            print_subsection(ki, omit=omit)
        for kj, vj in vi.items():
            if 'part' not in kj:
                print_ncol_loop(kj, vj, r, R, s, show=show, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)
                if show == 'u_R' or show == 'u_R_short':
                    print(omit)
                    return
            else:
                print_subsection(kj, omit=omit)
                for km, vm in vj.items():
                    print_ncol_loop(km, vm, r, R, s, show=show, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)

    print_section('q  space', omit=omit)
    for ki, vi in u['func_q'].items():
        print_subsection(ki, omit=omit)
        for kj, vj in vi.items():
            if 'part' not in kj:
                print_ncol(q, vj, header=['q', kj], show=show, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)
            else:
                print_subsection(kj, omit=omit)
                for km, vm in vj.items():
                    print_ncol(q, vm, header=['q', km], show=show, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)

    print(omit)

def print_func(u, r, name='unnamed', header=None, show='u_R_short', info=True, ncol=4,
               fmt=None, omit='', data_format=0):

    header = ['r', 'f(r)'] if header == None else header
    if info:
        print(TABLE_HEADER)
        print( TABLE_UNDERLINE)
        print_func_info(name, u.info)
    print_ncol(r, u(), header=header, show=show, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)

def print_all(u, r, q, R=None, s=None, title='', show='u_R_short', info=True, file=None,
              logo=False, ncol=4, fmt=None, omit='', data_format=0):
    '''
    It merges print_bifold(), print_bifold_di() and print_bifold_dd().
    And, if a file name given to "file" variable all output will be written to the file.

    show == 'info' only show names, norm, renorm, vol2, vol4, msr
    show == 'u_R' show the folded potential v.s. r including 'info'
    show == 'u_R_short' show the folded potential (summarized) v.s. r including 'info'
    show == 'all' show everything calculated
    show == 'all_short' show everything calculated in a summarized format
    '''

    if file is not None:
        import sys
        orig_stdout = sys.stdout
        f = open(file, 'w', encoding="utf-8")
        sys.stdout = f

    if logo:
        print_bifold_logo(omit=omit)

    if title!='':
        print_title(title=title, omit=omit)

    QUIT = False
    try:
        name = u['func_i']['u_R'][0]['name']
    except KeyError:
        name = u['func_i']['total']['u_R'][0]['name']
    except TypeError:
        name = 'unnamed'

    if name in ['u_direct', 'u_exchange_zr', 'u_coul_bifold_d']:
        print_bifold(u, r, q, R=R, s=s, show=show, info=info, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)
    elif name in ['u_bifold_zr', 'u_m3y_reid_zr', 'u_m3y_paris_zr']:
        print_bifold_di(u, r, q, R=R, s=s, show=show, info=info, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)
    elif 'dim3y' in name or 'ddm3y' in name or 'bdm3y' in name or 'cdm3y' in name:
        print_bifold_dd(u, r, q, R=R, s=s, show=show, info=info, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)
    elif name == 'unnamed':
        print_func(u, r, show=show, info=info, ncol=ncol, fmt=fmt, omit=omit, data_format=data_format)
    else:
        print(f'Something went wrong while printing! {name} is undefined.')
        QUIT = True

    if file is not None:
        sys.stdout = orig_stdout
        f.close()

    if QUIT:
        quit()


def print_potentials(u, r, q, R=None, s=None, title='', show='u_R', info=False, file=None,
                  logo=False, ncol=1, fmt=None, omit='#', data_format=0):

    fmt = ['>9.4f', '>17.9e'] if fmt == None else fmt

    print_all(u, r, q, R=R, s=s, title=title, show=show, info=info, file=file,
              logo=logo, ncol=ncol, fmt=fmt, omit=omit,
              data_format=data_format)