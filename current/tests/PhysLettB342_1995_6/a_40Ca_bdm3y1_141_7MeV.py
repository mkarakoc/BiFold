from bifold import *

r = mesh(zero, 10, 0.05)  # fm
q = mesh(zero,  3, 0.02)  # fm^-1
R = r.copy()
s = r.copy()

e_lab = 141.7
z_proj, a_proj =  2,  4
z_targ, a_targ = 20, 40

rho_p = f_2prm_gaussian(r, 0.4229, (1/0.7024)**.5)
rho_t = f_2prm_fermi(r, 0.169, 3.60, 0.523)

rc = 1.4 * (power(a_proj, 1 / 3) + power(a_targ, 1 / 3))
u_coul = u_coul_ucs(R, rc, z_proj, z_targ)

u_paris = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                  r, q, R, s, Cs=1/4, dd_name='bdm3y1', vnn_name='paris')

u_reid = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                  r, q, R, s, Cs=1/4, dd_name='bdm3y1', vnn_name='reid')

# Physics Letters B 342(1995) 6 - 12
# D.T. Khoa, W. von Oertzen.
# Fig 1. in page 8.
u_fig1_total_poly = [0.000949921, 0.00638314, -0.336656, 1.92118, 0.617326, 1.09765, -109.966]
u_fig1_paris_d_poly = [0.00532896, -0.143045, 1.44447, -6.64212, 13.3265, -11.0381, 3.94755, 52.221]
u_fig1_paris_exch_poly = [-0.00686376, 0.187661, -1.92189, 8.92557, -18.5332, 20.6983, -7.33647, -162.125]
u_fig1_reid_d_poly = [0.00253481, -0.0669833, 0.67631, -3.19149, 6.54708, -2.6208, 2.35148, -46.1083]
u_fig1_reid_exch_poly = [-0.00336784, 0.0926127, -0.958379, 4.53596, -9.6957, 10.5354, -4.88095, -62.8037]

u_fig1_total = f_extrapolate(R, 0.02, 7.7863, u_fig1_total_poly)
u_fig1_paris_d = f_extrapolate(R, 0.015, 7.609, u_fig1_paris_d_poly)
u_fig1_paris_exch = f_extrapolate(R, 0.0257, 7.6994, u_fig1_paris_exch_poly)
u_fig1_reid_d = f_extrapolate(R, 0.072, 7.5989, u_fig1_reid_d_poly)
u_fig1_reid_exch = f_extrapolate(R, 0.0961, 7.6895, u_fig1_reid_exch_poly)


print_all(u_paris, r, q)
print_all(u_reid, r, q)

markevery = int(len(R)*2.5/100)
figures = [
pplot(R, u_fig1_paris_d, label='Fig1-paris-direct', linestyle = 'None', marker='s', markevery=markevery, color='green', alpha=0.5),
pplot(R, u_fig1_reid_d, label='Fig1-reid-direct', linestyle = 'None', marker='v', markevery=markevery, color='blue', alpha=0.5),
pplot(R, u_fig1_paris_exch, label='Fig1-paris-exchange', linestyle = 'None', marker='^', markevery=markevery, color='red', alpha=0.5),
pplot(R, u_fig1_reid_exch, label='Fig1-reid-exchange', linestyle = 'None', marker='d', markevery=markevery, color='magenta', alpha=0.5),
pplot(R, u_fig1_total, label='Fig1-total', linestyle = 'None', marker='*', markevery=markevery, color='black', alpha=0.5)
    ]

plot_potentials([u_paris, u_reid], R, part='all', add_plot=figures)
