from bifold import *

r = mesh(zero, 10, 0.05)  # fm
q = mesh(zero,  3, 0.05)  # fm^-1

e_lab = 141.7 # MeV
a_proj =  4

rho_p = f_2prm_gaussian(r, 0.4229, (1/0.7024)**.5)
rho_t = f_2prm_fermi(r, 0.169, 3.60, 0.523)

u1 =  u_m3y_reid_zr(e_lab, a_proj, rho_p, rho_t, r, q)

print_all(u1, r, q)
plot_potentials(u1, r, part='all')
