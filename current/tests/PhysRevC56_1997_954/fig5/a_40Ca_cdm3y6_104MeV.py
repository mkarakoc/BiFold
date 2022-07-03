from bifold import *

r_min, r_max, dr = 5e-2, 15, 5/100 #fm
q_min, q_max, dq = 5e-2,  3, 5/100 #fm^-1

#r_min, r_max, dr = 4/100, 15, 4/100 #fm
#q_min, q_max, dq = 4/100,  3, 4/100 #fm^-1


r = mesh(r_min, r_max, dr)  # fm
q = mesh(q_min, q_max, dq)  # fm^-1
R = r.copy()
s = r.copy()

# Phys. Rev. C 56 (1997) 954.
# Dao T. Khoa, G. R. Satchler, and W. von Oertzen
u_fig5_cdm3y6_t_poly = [0.00259902, -0.0374923, 0.174945, -0.38456, 0.673565, 3.4871, 0.0990434, -122.409]
u_fig5_cdm3y4_t_poly = [0.00325549, -0.048278, 0.258383, -0.772158, 1.67862, 2.46336, 0.839039, -127.58]
u_fig5_cdm3y2_t_poly = [0.0115437, -0.168472, 0.724117, -1.06224, 5.24953, -0.221292, -131.727]
u_fig5_ddm3y1_t_poly = [0.00922679, -0.14986, 0.706704, -1.26995, 6.15651, -0.737525, -136.596]

u_fig5_cdm3y6_t = f_extrapolate(R, 0.01651, 5.0512, u_fig5_cdm3y6_t_poly)
u_fig5_cdm3y4_t = f_extrapolate(R, 0.10307, 4.98569, u_fig5_cdm3y4_t_poly)
u_fig5_cdm3y2_t = f_extrapolate(R, 0.10339, 5.0490, u_fig5_cdm3y2_t_poly)
u_fig5_ddm3y1_t = f_extrapolate(R, 0.13798, 5.02917, u_fig5_ddm3y1_t_poly)

# BiFold calculations
e_lab = 104.0
z_proj, a_proj =  2,  4
z_targ, a_targ = 20, 40
Cs = '1/4'

rho_p = f_2prm_gaussian(r, 0.4229, (1/0.7024)**.5)
rho_t = f_2prm_fermi(r, 0.169, 3.60, 0.523)

rc = 1.405 * (power(a_proj, 1 / 3) + power(a_targ, 1 / 3))
u_coul = u_coul_ucs(R, rc, z_proj, z_targ)

u_cdm3y6 = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                  r, q, R, s, Cs=eval(Cs), dd_name='cdm3y6', vnn_name='paris', u_ex_iter=8)

u_cdm3y4 = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                  r, q, R, s, Cs=eval(Cs), dd_name='cdm3y4', vnn_name='paris', u_ex_iter=8)

u_cdm3y2 = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                  r, q, R, s, Cs=eval(Cs), dd_name='cdm3y2', vnn_name='paris', u_ex_iter=8)

u_ddm3y1 = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                  r, q, R, s, Cs=eval(Cs), dd_name='ddm3y1', vnn_name='paris', u_ex_iter=8)

print_all(u_ddm3y1, r, q)
print_all(u_cdm3y2, r, q)
print_all(u_cdm3y4, r, q)
print_all(u_cdm3y6, r, q)

markevery = int(len(R)*2.5/100)
figures = [
pplot(R, u_fig5_cdm3y6_t, label='total - Fig 5. - 6',linestyle = 'None', marker='s', markevery=markevery, color='green', alpha=0.5),
pplot(R, u_fig5_cdm3y4_t, label='total - Fig 5. - 4',linestyle = 'None', marker='v', markevery=markevery, color='red', alpha=0.5),
pplot(R, u_fig5_cdm3y2_t, label='total - Fig 5. - 2',linestyle = 'None', marker='^', markevery=markevery, color='blue', alpha=0.5),
pplot(R, u_fig5_ddm3y1_t, label='total - Fig 5. - 1',linestyle = 'None', marker='d', markevery=markevery, color='blue', alpha=0.5)
    ]
suptitle = f'BiFold vs literature for CDM3Y6/Paris, Cs = {Cs}'
title = 'Reference: Fig. 5. @ Phys. Rev. C 56 (1997) 954.'

plot_potentials([u_cdm3y6, u_cdm3y4, u_cdm3y2, u_ddm3y1], R, part='all',
                add_plot=figures, xlimit=(0, 5), suptitle=suptitle, title=title)
