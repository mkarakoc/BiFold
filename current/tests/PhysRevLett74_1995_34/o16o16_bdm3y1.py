from bifold import *

r = mesh(zero, 10, 0.05)  # fm
q = mesh(zero, 3, 0.05)  # fm^-1
R = r.copy()  # fm
s = r.copy()

# Phys. Rev. Lett. 74 (1995) 34.
# Dao T. Khoa, W. von Oertzen, H. G. Bohlen, G. Bartnitzky, H. Clement, Y. Sugiyama,
# B. Gebauer, A. N. Ostrowski, Th. Wilpert, M. Wilpert, and C. Langner
u_fig1_bdm3y1_d_poly = [0.0034318, -0.089972, 0.915882, -4.90944, 17.7298, -44.6933, 5.53461, 223.865]
u_fig1_bdm3y1_e_poly = [-0.00332156, 0.0859079, -0.831079, 4.19002, -17.706, 65.3812, -6.13047, -553.627]
u_fig1_bdm3y1_t_poly = [-0.00280904, 0.0734505, -0.704694, 3.15484, -9.48535, 31.7323, -4.97658, -330.222]

u_fig1_bdm3y1_d = f_extrapolate(R, 0.0253, 7.9544, u_fig1_bdm3y1_d_poly)
u_fig1_bdm3y1_e = f_extrapolate(R, 0.0122, 7.9904, u_fig1_bdm3y1_e_poly)
u_fig1_bdm3y1_t = f_extrapolate(R, 0.0101, 7.9896, u_fig1_bdm3y1_t_poly)

e_lab = 250.0
z_proj, a_proj =  8, 16
z_targ, a_targ =  8, 16

rho_p = f_2prm_fermi(r, 0.181, 2.525, 0.450)
rho_t = f_2prm_fermi(r, 0.181, 2.525, 0.450)

rc = 1.405 * (power(a_proj, 1 / 3) + power(a_targ, 1 / 3))

u_coul = u_coul_ucs(R, rc, z_proj, z_targ)

u_bdm3y1 = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                  r, q, R, s, Cs=1/4, dd_name='bdm3y1', vnn_name='paris')

print_all(u_bdm3y1, r, q)

markevery = int(len(R)*2.5/100)
fig1 = [
pplot(R, u_fig1_bdm3y1_d, label='direct - Fig 1.',linestyle = 'None', marker='^', markevery=markevery, color='red', alpha=0.5),
pplot(R, u_fig1_bdm3y1_e, label='exchange - Fig 1.',linestyle = 'None', marker='v', markevery=markevery, color='blue', alpha=0.5),
pplot(R, u_fig1_bdm3y1_t, label='total - Fig 1.',linestyle = 'None', marker='s', markevery=markevery, color='green', alpha=0.5)
    ]

suptitle = 'BiFold vs literature for BDM3Y1/Paris'
title = 'Reference: Fig. 1. @ Phys. Rev. Lett. 74 (1995) 34.'
plot_potentials(u_bdm3y1, R, part='all', add_plot=fig1, suptitle=suptitle, title=title,
                linestyles=['dashed', 'dotted', 'solid'], xlimit=(0,8))
