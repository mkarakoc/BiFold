from bifold import *

e_lab = 141.7
z_proj, a_proj =  2,  4
z_targ, a_targ = 20, 40

r = mesh(zero, 10, 0.05)  # fm
q = mesh(zero,  3, 0.05)  # fm^-1

rho_p = f_2prm_gaussian(r, 0.4229, (1/0.7024)**.5)
rho_t = f_2prm_fermi(r, 0.169, 3.60, 0.523)
#rho_t = f_nudat(r, z_targ, a_targ)
u = u_ddm3y_reid_zr(e_lab, a_proj, rho_p, rho_t, r, q)

print_all(u, r, q)
plot_potentials(u, r, part='all')
