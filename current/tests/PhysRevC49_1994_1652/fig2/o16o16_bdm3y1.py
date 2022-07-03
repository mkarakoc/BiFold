from bifold import *

r = mesh(zero, 10, 0.05)  # fm
q = mesh(zero,  3, 0.02)  # fm^-1
R = r.copy()
s = r.copy()

# Phys. Rev. C 49 (1994) 1652.
u_fig1_bdm3y1_t_poly = [-0.00234749, 0.0614265, -0.581974, 2.52119, -7.61662, 28.3349, -4.34508, -317.04]
u_fig1_bdm3y1_t = f_extrapolate(R, 0.00063, 7.86524, u_fig1_bdm3y1_t_poly)

e_lab = 250.0
z_proj, a_proj =  8, 16
z_targ, a_targ =  8, 16

rho_p = f_2prm_fermi(r, 0.181, 2.525, 0.450)
rho_t = f_2prm_fermi(r, 0.181, 2.525, 0.450)

rc = 1.405 * (power(a_proj, 1 / 3) + power(a_targ, 1 / 3))
rc = 2.71*2
u_coul = u_coul_ucs(R, rc, z_proj, z_targ)

u_bdm3y1 = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                  r, q, R, s, Cs=1/4, dd_name='bdm3y1', vnn_name='reid')
print_all(u_bdm3y1, r, q)

markevery = int(len(R)*2.5/100)
fig2 = pplot(R, u_fig1_bdm3y1_t, label='total - Fig 2.',linestyle = 'None', marker='s',
             markevery=markevery, color='green', alpha=0.5)
suptitle = 'BiFold vs literature for BDM3Y1/Reid'
title = 'Reference: Fig. 2. @ Phys. Rev. C 49 (1994) 1652.'
plot_potentials(u_bdm3y1, R, part='all', add_plot=fig2, suptitle=suptitle, title=title,
                linestyles=['dashed', 'dotted', 'solid'])
