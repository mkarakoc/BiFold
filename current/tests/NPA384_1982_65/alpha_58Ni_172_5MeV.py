from bifold.filon import *

e_lab = 172.5
z_proj, a_proj =  2,  4
z_targ, a_targ = 28, 58

r = mesh(zero, 11, 0.02)  # fm
q = mesh(zero,  3, 0.02)  # fm^-1

rho_p = f_2prm_gaussian(r, 0.4229, (1/0.7024)**.5)
rho_t = f_nudat(r, z_targ, a_targ)
u_ddm3y = u_ddm3y_reid_zr(e_lab, a_proj, rho_p, rho_t, r, q)
u_m3y = u_m3y_reid_zr(e_lab, a_proj, rho_p, rho_t, r, q)

title = f'alpha + 58Ni@ {e_lab} NPA 384 (1982) 65'

print_all(u_ddm3y, r, q, title=f'DDM3Y - ZR - {title}')
print_all(u_m3y, r, q, title=f'M3Y - ZR - {title}')
#plot_potentials(u, r, part='all', semilogy=True)
plot_potentials([u_m3y, u_ddm3y], r, semilogy=True, ylimit=(0.01, 999), figsize=(6,7.1), title=title)
