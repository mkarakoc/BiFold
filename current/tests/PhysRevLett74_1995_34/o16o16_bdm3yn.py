from bifold.simpson import *

r = mesh(zero, 10, 0.05)  # fm
q = mesh(zero,  3, 0.05)  # fm^-1


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

u_fig1_bdm3y2_d_poly = [0.0071336, -0.191227, 1.98369, -10.2205, 29.5158, -52.8523, 15.2599, 182.156]
u_fig1_bdm3y2_t_poly = [-0.00205206, 0.0546904, -0.531231, 2.38045, -7.33398, 26.8773, -5.54677, -299.266]
u_fig1_bdm3y2_d = f_extrapolate(R, 0.0290, 7.8645, u_fig1_bdm3y2_d_poly)
u_fig1_bdm3y2_t = f_extrapolate(R, 0.0283, 7.9486, u_fig1_bdm3y2_t_poly)


u_fig1_bdm3y3_d_poly = [-0.00095538,0.034498,-0.512948,4.0515,-18.3576,48.487,-72.6893,50.8867,-14.6314,118.462]
u_fig1_bdm3y3_t_poly = [-0.00381729, 0.105935, -1.10422, 5.35587, -13.6016, 27.3881, -7.71518, -248.507]
u_fig1_bdm3y3_d = f_extrapolate(R, 0.0349, 7.6810, u_fig1_bdm3y3_d_poly)
u_fig1_bdm3y3_t = f_extrapolate(R, 0.0207, 7.9551, u_fig1_bdm3y3_t_poly)


e_lab = 250.0
z_proj, a_proj =  8, 16
z_targ, a_targ =  8, 16

rho_p = f_2prm_fermi(r, 0.181, 2.525, 0.450)
rho_t = f_2prm_fermi(r, 0.181, 2.525, 0.450)

rc = 1.405 * (power(a_proj, 1 / 3) + power(a_targ, 1 / 3))

u_coul = u_coul_ucs(R, rc, z_proj, z_targ)

u_bdm3y1 = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                  r, q, R, s, Cs=1/4, dd_name='bdm3y1', vnn_name='paris')

u_bdm3y2 = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                  r, q, R, s, Cs=1/4, dd_name='bdm3y2', vnn_name='paris')

u_bdm3y3 = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                  r, q, R, s, Cs=1/4, dd_name='bdm3y3', vnn_name='paris')

#u_bdm3y1 = u_xdm3yn_zr(e_lab, a_proj, rho_p, rho_t, r, q, R, s, dd_name='bdm3y1', vnn_name='paris')
#u_bdm3y2 = u_xdm3yn_zr(e_lab, a_proj, rho_p, rho_t, r, q, R, s, dd_name='bdm3y2', vnn_name='paris')
#u_bdm3y3 = u_xdm3yn_zr(e_lab, a_proj, rho_p, rho_t, r, q, R, s, dd_name='bdm3y3', vnn_name='paris')

print_all(u_bdm3y1, r, q)
print_all(u_bdm3y2, r, q)
print_all(u_bdm3y3, r, q)

markevery = int(len(R)*2.5/100)
fig1_d = [
pplot(R, u_fig1_bdm3y1_d, label='direct - bdm3y1',linestyle = 'None', marker='^', markevery=markevery, color='green', alpha=0.5, zorder=999),
pplot(R, u_fig1_bdm3y2_d, label='direct - bdm3y2',linestyle = 'None', marker='^', markevery=markevery, color='blue', alpha=0.5, zorder=999),
pplot(R, u_fig1_bdm3y3_d, label='direct - bdm3y3',linestyle = 'None', marker='^', markevery=markevery, color='red', alpha=0.5, zorder=999)
    ]

fig1_t = [
#pplot(R, u_fig1_bdm3y1_e, label='exchange - Fig 1.',linestyle = 'None', marker='v', markevery=markevery, color='blue', alpha=0.5, zorder=999),
pplot(R, u_fig1_bdm3y1_t, label='total - bdm3y1',linestyle = 'None', marker='s', markevery=markevery, color='green', alpha=0.5, zorder=999),
pplot(R, u_fig1_bdm3y2_t, label='total - bdm3y2',linestyle = 'None', marker='s', markevery=markevery, color='blue', alpha=0.5, zorder=999),
pplot(R, u_fig1_bdm3y3_t, label='total - bdm3y3',linestyle = 'None', marker='s', markevery=markevery, color='red', alpha=0.5, zorder=999)
    ]


suptitle = 'BiFold vs literature for BDM3Y1/Paris'
title = 'Reference: Fig. 1. @ Phys. Rev. Lett. 74 (1995) 34.'
plot_potentials([u_bdm3y1, u_bdm3y2, u_bdm3y3], R, part='direct', add_plot=fig1_d, suptitle=suptitle, title=title,
                linestyles=['dashed', 'dotted', 'solid'], xlimit=(0,8))

plot_potentials([u_bdm3y1, u_bdm3y2, u_bdm3y3], R, part='total', add_plot=fig1_t, suptitle=suptitle, title=title,
                linestyles=['dashed', 'dotted', 'solid'], xlimit=(0,8))


#plot_potentials([u_bdm3y1, u_bdm3y2, u_bdm3y3], R, part='direct', add_plot=fig1, suptitle=suptitle, title=title,
#                linestyles=['dashed', 'dotted', 'solid'], xlimit=(0,8))
