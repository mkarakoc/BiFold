from bifold import *
#print_bifold_logo()

r = mesh(zero, 10, 0.05)  # fm
q = mesh(zero,  3, 0.05)  # fm^-1
R = r.copy()  # fm
s = r.copy()

e_lab = 141.7
z_proj, a_proj =  2,  4
z_targ, a_targ = 20, 40

rho_p = f_2prm_gaussian(r, 0.4229, (1/0.7024)**.5)
#rho_t = f_2prm_fermi(r, 0.169, 3.60, 0.523)
rho_t = f_ripl(r, z_targ, a_targ)

u1 =  u_m3y_reid_zr(e_lab, a_proj, rho_p, rho_t, r, q, R, s)
u2 =  u_m3y_paris_zr(e_lab, a_proj, rho_p, rho_t, r, q, R, s)

print_all(u1, r, q)
print_all(u2, r, q, show='info')
plot_potentials([u1, u2], R, part='all', linestyles=['dashed', 'dotted', 'solid', 'dashed', 'dotted', 'solid'])
plot_fouriers([u1, u2], q, part='all', linestyles=['dashed', 'dotted', 'solid', 'dashed', 'dotted', 'solid'])
plot_all(u1, r, q)
