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
This module contains the plot tools for BiFold.
"""

from matplotlib import pyplot as plt
from random import randint

LINE_STYLES = ['solid', 'dashed', 'dashdot', 'dotted']
COLORS = 'krbgcmy'

def colors():
    return COLORS[randint(0, len(COLORS) - 1)]

def linestyles():
    return LINE_STYLES[randint(0, len(LINE_STYLES) - 1)]

def pplot(*args, **kwargs):
    return args, kwargs

def plot_potentials(u, R, part='total', space='func_r', func='u_R', legends=None, legend_ext=None, figsize=(9,6),
                    linestyles=None, colors=None, linewidths=None, loc='best', label_x='R [fm]',
                    label_y='U(R) [MeV]', file=None, show=True, block=True, title=False, suptitle='',
                    semilogy=False, xlimit=None, ylimit=None, add_plot=None):

    # to plot all the given potentials with the all parts
    if part == 'all' and isinstance(u, list):
        if legend_ext[0] != '':
            legend_ext =  [[li, li, li] for li in legend_ext]
            legend_ext = sum(legend_ext, [])

        part = [['direct', 'exchange', 'total'] for ui in u]
        part = sum(part, [])
        u = [[ui, ui, ui] for ui in u]
        u = sum(u, [])

    # to plot all the given potential with the all parts
    u = [u,u,u] if part == 'all' else u
    part = ['direct', 'exchange', 'total'] if part == 'all' else part

    # to plot only one part of a potential
    u = u if isinstance(u, list) else [u]
    part = part if isinstance(part, list) else [part]
    part = part if len(u)==len(part) else [part[0]] * len(u)

    linestyles =  ['solid', 'dashed', 'dashdot', 'dotted'] if linestyles==None else linestyles
    colors = 'krbgcmy' if colors==None else colors
    linewidths = [2]*len(u) if linewidths==None else linewidths

    plt.figure(figsize=figsize)

    if suptitle!='':
        plt.suptitle(suptitle)

    if isinstance(title, str):
        plt.title(title)
    elif title:
        if len(part)>1:
            part_str = ' - '.join([part_i for part_i in part])
        else:
            part_str = part[0]
        plt.title(part_str)

    if add_plot != None:
        add_plot = add_plot if isinstance(add_plot, list) else [add_plot]
        for add_i in add_plot:
            plt.plot(*add_i[0], **add_i[1])

    plt.grid(color='lightgray', zorder=-999, linestyle='dashed')
    for ni, ui in enumerate(u):
        try:
            u_R = ui[space][func]
        except KeyError:
            u_R = ui[space][part[ni]][func]
        except TypeError:
            u_R = ui()

        if legends==None:
            try:
                name = ui['func_i'][func][0]['name']
            except KeyError:
                name = ui['func_i'][part[ni]][func][0]['name']
                if legend_ext!=None:
                     name += legend_ext[ni]
            except TypeError:
                name = ui.info[0]['name']
        else:
            name = legends[ni]

        if semilogy:
            u_R = abs(u_R)
            plt.semilogy(R, u_R, label=name, linestyle=linestyles[ni%len(linestyles)],
                     color=colors[ni%len(colors)], linewidth=linewidths[ni%len(linewidths)])
            plt.ylabel(f'|{label_y}|', fontsize=14)
        else:
            plt.plot(R, u_R, label=name, linestyle=linestyles[ni%len(linestyles)],
                     color=colors[ni%len(colors)], linewidth=linewidths[ni%len(linewidths)])
            plt.ylabel(label_y, fontsize=14)


    plt.legend(loc=loc)
    plt.xlabel(label_x, fontsize=14)


    plt.tick_params(direction='in', labelsize=14, length=7, width=1)
    plt.minorticks_on()
    plt.tick_params(which='minor', direction='in', length=5, width=1)

    if xlimit==None:
        plt.xlim(R[0], R[-1])
    else:
        plt.xlim(*xlimit)

    if ylimit!=None:
        plt.ylim(*ylimit)

    if file!=None:
        if file[0]=='.':
            plt.savefig(f'{name}{file}')
        else:
            plt.savefig(file)

    if show:
        plt.show(block=block)



def plot_fouriers(u, q, part='total', space='func_q', func='u_R', legends=None, legend_ext=None, figsize=(9,6),
                  linestyles=None, colors=None, linewidths=None, loc='best',
                  label_x='q [fm$^{-1}$]', label_y='U(q) [MeV fm$^3$]', file=None, show=True, block=True, title=False,
                  suptitle='', semilogy=False, xlimit=None, ylimit=None, add_plot=None):
    plot_potentials(u, q, part=part, space=space, func=func, legends=legends, legend_ext=legend_ext, figsize=figsize,
                    linestyles=linestyles, colors=colors, linewidths=linewidths, loc=loc,
                    label_x=label_x, label_y=label_y, file=file, show=show, block=block, title=title,
                    suptitle=suptitle, semilogy=semilogy, xlimit=xlimit, ylimit=ylimit, add_plot=add_plot)


def plot_bifold(u, r, q, R=None, s=None):
    '''
    It plots only direct/exchange part of the density independent folded potentials.
    '''

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    fig, ax = plt.subplots(3, 2, figsize=(12, 9))

    for k, v in u['func_r'].items():
        if k == 'u_R':
            label = u['func_i']['u_R'][0]['name']
            ax[0,0].plot(R, v, label=label, color=colors(), linestyle=linestyles())
            ax[0,0].set_xlim(R[0], R[-1])
        elif 'rho_' in k:
            ax[1,0].plot(r, v, label=k, color=colors(), linestyle=linestyles())
            ax[1,0].set_xlim(r[0], r[-1] / 2)
        elif 'vnn' in k:
            ax[2,0].text(.5, .9, '$v_{nn}$', horizontalalignment='center', transform=ax[2,0].transAxes)
            ax[2,0].plot(s, v, label=k, color=colors(), linestyle=linestyles())
            ax[2,0].set_xlim(s[0], s[-1] / 2)
            v_min = min(v)
            ax[2,0].set_ylim(v_min * 1.1, abs(v_min) * 2)

    for k, v in u['func_q'].items():
        if k == 'u_R':
            label = u['func_i']['u_R'][0]['name']
            ax[0,1].plot(q, v, label=label, color=colors(), linestyle=linestyles())
            ax[0,1].set_xlim(q[0], q[-1])
        elif 'rho_' in k:
            ax[1,1].plot(q, v, label=k, color=colors(), linestyle=linestyles())
            ax[1,1].set_xlim(q[0], q[-1])
        elif 'vnn' in k:
            ax[2,1].text(.5, .9, '$v_{nn}$', horizontalalignment='center', transform=ax[2,1].transAxes)
            ax[2,1].plot(q, v, label=k, color=colors(), linestyle=linestyles())
            ax[2,1].set_xlim(q[0], q[-1])
            v_min = min(v)
            ax[2,1].set_ylim(v_min * 1.1, abs(v_min) * 2)


    for i in range(3):
        for j in range(2):
            ax[i,j].grid(color='lightgray', zorder=-999, linestyle='dashed')
            if not ax[i,j].lines:
                ax[i,j].set_visible(False)
            else:
                ax[i,j].legend(loc='best')
                ax[i,j].xaxis.set_ticks_position('none')
                ax[i,j].yaxis.set_ticks_position('none')

    plt.show()



def plot_bifold_di(u, r, q, R=None, s=None):
    '''
    It plots direct and exchange parts of the density independent folded potentials.
    '''

    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    ax[0,1].text(.5,1.1,'direct',
        horizontalalignment='center',
        transform=ax[0,1].transAxes)

    ax[0,2].text(.5,1.1,'exchange',
        horizontalalignment='center',
        transform=ax[0,2].transAxes)

    for ki, vi in u['func_r'].items():
        for kj, vj in vi.items():
            if kj == 'u_R':
                label = u['func_i'][ki]['u_R'][0]['name']
                ax[0, 0].plot(R, vj, label=label, color=colors(), linestyle=linestyles())
                ax[0, 0].set_xlim(R[0], R[-1])
            elif  'rho_'in kj and 'direct' in ki:
                ax[0,1].plot(r, vj, label=kj, color=colors(), linestyle=linestyles())
                ax[0, 1].set_xlim(r[0], r[-1] / 2)
            elif  'rho_'in kj and 'exchange' in ki:
                ax[0,2].plot(r, vj, label=kj, color=colors(), linestyle=linestyles())
                ax[0, 2].set_xlim(r[0], r[-1] / 2)
            elif 'vnn' in kj and 'direct' in ki:
                ax[1,1].text(.5, .9, '$v_{nn}$', horizontalalignment='center',
                              transform=ax[1, 1].transAxes)
                ax[1,1].plot(s, vj, label=ki, color=colors(), linestyle=linestyles())
                ax[1, 1].set_xlim(s[0], s[-1] / 2)
                vj_min = min(vj)
                ax[1,1].set_ylim(vj_min * 1.1, abs(vj_min) * 2)
            elif 'vnn' in kj and 'exchange' in ki:
                ax[1,2].text(.5, .9, '$v_{nn}$', horizontalalignment='center',
                              transform=ax[1, 2].transAxes)
                ax[1,2].plot(s, vj, label=ki, color=colors(), linestyle=linestyles())
                ax[1, 2].set_xlim(s[0], s[-1] / 2)
                vj_min = min(vj)
                ax[1,2].set_ylim(vj_min * 1.1, abs(vj_min) * 2)


    for i in range(2):
        for j in range(3):
            ax[i,j].grid(color='lightgray', zorder=-999, linestyle='dashed')
            if not ax[i,j].lines:
                ax[i,j].set_visible(False)
            elif i==3 and j==0:
                ax[i,j].legend(loc='upper right')
                ax[i,j].xaxis.set_ticks_position('none')
                ax[i,j].yaxis.set_ticks_position('none')
            else:
                ax[i,j].legend(loc='best')
                ax[i,j].xaxis.set_ticks_position('none')
                ax[i,j].yaxis.set_ticks_position('none')

    plt.show()

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    ax[0,1].text(.5,1.1,'direct',
        horizontalalignment='center',
        transform=ax[0,1].transAxes)

    ax[0,2].text(.5,1.1,'exchange',
        horizontalalignment='center',
        transform=ax[0,2].transAxes)

    for ki, vi in u['func_q'].items():
        for kj, vj in vi.items():
            if kj == 'u_R':
                label = u['func_i'][ki]['u_R'][0]['name']
                ax[0, 0].plot(q, vj, label=label, color=colors(), linestyle=linestyles())
            elif  'rho_'in kj and 'direct' in ki:
                ax[0,1].plot(q, vj, label=kj, color=colors(), linestyle=linestyles())
            elif  'rho_'in kj and 'exchange' in ki:
                ax[0,2].plot(q, vj, label=kj, color=colors(), linestyle=linestyles())
            elif  'vnn'in kj and 'direct' in ki:
                ax[1,1].text(.5, .9, '$v_{nn}$', horizontalalignment='center',
                              transform=ax[1, 1].transAxes)
                ax[1,1].plot(q, vj, label=ki, color=colors(), linestyle=linestyles())
            elif  'vnn'in kj and 'exchange' in ki:
                ax[1,2].text(.5, .9, '$v_{nn}$', horizontalalignment='center',
                              transform=ax[1, 2].transAxes)
                ax[1,2].plot(q, vj, label=ki, color=colors(), linestyle=linestyles())

    for i in range(2):
        for j in range(3):
            ax[i,j].grid(color='lightgray', zorder=-999, linestyle='dashed')
            if not ax[i,j].lines:
                ax[i,j].set_visible(False)
            elif i==3 and j==0:
                ax[i,j].legend(loc='upper right')
                ax[i,j].xaxis.set_ticks_position('none')
                ax[i,j].yaxis.set_ticks_position('none')
                ax[i, j].set_xlim(q[0], q[-1])
            else:
                ax[i,j].legend(loc='best')
                ax[i,j].xaxis.set_ticks_position('none')
                ax[i,j].yaxis.set_ticks_position('none')
                ax[i, j].set_xlim(q[0], q[-1])
    plt.show()



def plot_bifold_dd(u, r, q, R=None, s=None):
    '''
    It plots direct and exchange parts of the density dependent folded potentials.
    '''


    R = r.copy() if R is None else R
    s = r.copy() if s is None else s

    fig, ax =  plt.subplots(4, 3, figsize=(16,9))
    ax[3, 0].text(.5, .9, '$v_{nn}$',
                  horizontalalignment='center',
                  transform=ax[3, 0].transAxes)

    ax[0,1].text(.5,1.1,'direct',
        horizontalalignment='center',
        transform=ax[0,1].transAxes)

    ax[0,2].text(.5,1.1,'exchange',
        horizontalalignment='center',
        transform=ax[0,2].transAxes)


    for ki, vi in u['func_r'].items():
        exchange_p = 1 if 'exchange' in ki else 0
        for kj, vj in vi.items():
            if 'part' not in kj:
                if kj=='u_R':
                    label = u['func_i'][ki]['u_R'][0]['name']
                    ax[0,0].plot(R, vj, label=label, color=colors(), linestyle=linestyles())
                    ax[0,0].set_xlim(R[0], R[-1])
                elif  'u_coul' in kj:
                    ax[2,0].plot(R, vj, label=kj, color=colors(), linestyle=linestyles())
                    ax[2,0].set_xlim(R[0], R[-1])
                elif  'rho_'in kj:
                    ax[0,2].plot(r, vj, label=kj, color=colors(), linestyle=linestyles())
                    ax[0,2].set_xlim(r[0], r[-1] / 2)
                elif 'vnn' in kj:
                    ax[1,2].text(.5, .9, '$v_{nn}$',
                                  horizontalalignment='center',
                                  transform=ax[1, 2].transAxes)
                    ax[1,2].plot(s, vj, label=ki, color=colors(), linestyle=linestyles())
                    ax[1,2].set_xlim(s[0], s[-1] / 2)
                    vj_min = min(vj)
                    ax[1,2].set_ylim(vj_min * 1.1, abs(vj_min) * 2)
                elif 'kf_' in kj:
                    ax[2,2].plot(r, vj, label=kj, color=colors(), linestyle=linestyles())
                    ax[2,2].set_xlim(r[0], r[-1])
            else:
                for km, vm in vj.items():
                    if km == 'u_R':
                        which_part = 'direct parts \n of $U_R$' if exchange_p==0 else 'exchange parts \n of $U_R$'
                        ax[1 + exchange_p, 0].text(.55, .5, which_part,
                                      horizontalalignment='center',
                                      transform=ax[1 + exchange_p, 0].transAxes)
                        ax[1 + exchange_p, 0].plot(R, vm, label=kj, color=colors(), linestyle=linestyles())
                        ax[1 + exchange_p, 0].set_xlim(R[0], R[-1])

                    if 'rho_' in km:
                        j = int(kj[-1]) - 1
                        ax[j, 1 + exchange_p].plot(r, vm, label=km, color=colors(), linestyle=linestyles())
                        ax[j, 1 + exchange_p].text(.5, .9, kj,
                                      horizontalalignment='center',
                                      transform=ax[j, 1 + exchange_p].transAxes)
                        ax[j, 1+ exchange_p].set_xlim(r[0], r[-1] / 2)

                    if 'vnn' in km:
                        j = int(kj[-1])
                        if j == 1:
                            ax[3, 0].plot(s, vm, label=ki, color=colors(), linestyle=linestyles())
                            ax[3, 0].set_xlim(s[0], s[-1]/2)
                            vm_min = min(vm)
                            ax[3, 0].set_ylim(vm_min*1.1, abs(vm_min)*2)


    for i in range(4):
        for j in range(3):
            ax[i,j].grid(color='lightgray', zorder=-999, linestyle='dashed')
            if not ax[i,j].lines:
                ax[i,j].set_visible(False)
            elif i==3 and j==0:
                ax[i,j].legend(loc='upper right')
                ax[i,j].xaxis.set_ticks_position('none')
                ax[i,j].yaxis.set_ticks_position('none')
            else:
                ax[i,j].legend(loc='best')
                ax[i,j].xaxis.set_ticks_position('none')
                ax[i,j].yaxis.set_ticks_position('none')

    plt.show()


    fig, ax =  plt.subplots(4, 3, figsize=(16,9))
    ax[3, 0].text(.5, .9, '$v_{nn}$',
                  horizontalalignment='center',
                  transform=ax[3, 0].transAxes)

    ax[0,1].text(.5,1.1,'direct',
        horizontalalignment='center',
        transform=ax[0,1].transAxes)

    ax[0,2].text(.5,1.1,'exchange',
        horizontalalignment='center',
        transform=ax[0,2].transAxes)

    for ki, vi in u['func_q'].items():
        exchange_p = 1 if 'exchange' in ki else 0
        for kj, vj in vi.items():
            if 'part' not in kj:
                if kj=='u_R':
                    label = u['func_i'][ki]['u_R'][0]['name']
                    ax[0,0].plot(q, vj, label=label, color=colors(), linestyle=linestyles())
                elif  'rho_'in kj:
                    ax[0,2].plot(q, vj, label=kj, color=colors(), linestyle=linestyles())
                elif  'vnn'in kj:
                    ax[1, 2].text(.5, .9, '$v_{nn}$',
                                  horizontalalignment='center',
                                  transform=ax[1, 2].transAxes)
                    ax[1,2].plot(q, vj, label=ki, color=colors(), linestyle=linestyles())
                elif  'kf_'in kj:
                    ax[2,2].plot(q, vj, label=kj, color=colors(), linestyle=linestyles())
            else:
                for km, vm in vj.items():
                    if km == 'u_R':
                        which_part = 'direct parts \n of $U_R$' if exchange_p==0 else 'exchange parts \n of $U_R$'
                        ax[1 + exchange_p, 0].text(.55, .5, which_part,
                                      horizontalalignment='center',
                                      transform=ax[1 + exchange_p, 0].transAxes)
                        ax[1 + exchange_p, 0].plot(q, vm, label=kj, color=colors(), linestyle=linestyles())

                    if 'rho_' in km:
                        j = int(kj[-1]) - 1
                        ax[j, 1 + exchange_p].plot(q, vm, label=km, color=colors(), linestyle=linestyles())
                        ax[j, 1 + exchange_p].text(.5, .9, kj,
                                      horizontalalignment='center',
                                      transform=ax[j, 1 + exchange_p].transAxes)
                    if 'vnn' in km:
                        j = int(kj[-1])
                        if j == 1:
                            ax[3, 0].plot(q, vm, label=ki, color=colors(), linestyle=linestyles())



    for i in range(4):
        for j in range(3):
            ax[i,j].grid(color='lightgray', zorder=-999, linestyle='dashed')
            if not ax[i,j].lines:
                ax[i,j].set_visible(False)
            elif i==3 and j==0:
                ax[i,j].legend(loc='upper right')
                ax[i,j].xaxis.set_ticks_position('none')
                ax[i,j].yaxis.set_ticks_position('none')
                ax[i,j].set_xlim(q[0], q[-1])
            else:
                ax[i,j].legend(loc='best')
                ax[i,j].xaxis.set_ticks_position('none')
                ax[i,j].yaxis.set_ticks_position('none')
                ax[i, j].set_xlim(q[0], q[-1])
    #fig.subplots_adjust(hspace=0.0)
    plt.show()



def plot_all(u, r, q, R=None, s=None, title=''):
    '''
    It merges plot_bifold(), plot_bifold_di() and plot_bifold_dd().
    And, if a file name given to "file" variable all output will be written to the file.
    '''

    if title!='':
        print_title(title=title, omit=omit)

    QUIT = False
    try:
        name = u['func_i']['u_R'][0]['name']
    except KeyError:
        name = u['func_i']['total']['u_R'][0]['name']

    if name in ['u_direct', 'u_exchange_zr']:
        plot_bifold(u, r, q, R=R, s=s)
    elif name in ['u_bifold_zr', 'u_m3y_reid_zr', 'u_m3y_paris_zr']:
        plot_bifold_di(u, r, q, R=R, s=s)
    elif 'dim3y' in name or 'ddm3y' in name or 'bdm3y' in name or 'cdm3y' in name:
        plot_bifold_dd(u, r, q, R=R, s=s)
    else:
        print(f'Something went wrong while printing! {name} is undefined.')
        QUIT = True

    if QUIT:
        quit()
