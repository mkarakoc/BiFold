from bifold import *

r_min, r_max, dr = zero, 15, 5/100 #fm
q_min, q_max, dq = zero,  3, 5/100 #fm^-1

r = mesh(r_min, r_max, dr)  # fm
q = mesh(q_min, q_max, dq)  # fm^-1
R = r.copy()
s = r.copy()

# Phys. Rev. C 63 (2001) 034007.
u_fig2_cdm3y6_t_poly = [0.000756501,-0.0289373,0.416498,-2.74808,7.49829,-3.00499,5.82575,-86.2282]
u_fig2_cdm3y6_t = f_extrapolate(R, 0.1295, 7.9639, u_fig2_cdm3y6_t_poly)

e_lab = 120.0
z_proj, a_proj =  2,  4
z_targ, a_targ =  6, 12
Cs = '1/4'

rho_p = f_2prm_gaussian(r, 4/(sqrt_pi*1.2658)**3, 1.2658)
rho_t = f_2prm_fermi(r, 0.19354, 2.214, 0.425)

rc = 1.405 * (power(a_proj, 1 / 3) + power(a_targ, 1 / 3))
u_coul = u_coul_ucs(R, rc, z_proj, z_targ)

u_cdm3y6 = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                  r, q, R, s, Cs=eval(Cs), dd_name='cdm3y6', vnn_name='paris')

print_all(u_cdm3y6, r, q)

markevery = int(len(R)*2.5/100)
fig2 = pplot(R, u_fig2_cdm3y6_t, label='total - Fig 5.',linestyle = 'None', marker='s', markevery=markevery, color='green', alpha=0.5)
suptitle = f'BiFold vs literature for CDM3Y6/Paris, Cs = {Cs}'
title = 'Reference: Fig. 2. @ Phys. Rev. C 63 (2001) 034007.'
plot_potentials(u_cdm3y6, R, part='all', add_plot=fig2, xlimit=(0, 8), suptitle=suptitle, title=title,
                linestyles=['dashed', 'dotted', 'solid'])
