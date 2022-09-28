from bifold import *

e_lab = 100.0
z_proj, a_proj =  2, 3
z_targ, a_targ =  2, 4

r = mesh(zero, 15, 0.02)  # fm
q = mesh(zero, 15, 0.02)  # fm^-1

rho_p = f_2prm_gaussian(r, 0.2201, (1/0.5505)**.5)
rho_t = f_2prm_gaussian(r, 0.4229, (1/0.7024)**.5)

u_ddm3y = u_ddm3y_reid_zr(e_lab, a_proj, rho_p, rho_t, r, q)
u_m3y = u_m3y_reid_zr(e_lab, a_proj, rho_p, rho_t, r, q)

title = f'3He + 4He @ {e_lab} MeV'
plot_potentials([u_m3y, u_ddm3y], r,  title=title)

plot_potentials([u_m3y, u_ddm3y], r,  title=title, semilogy=True)

print_all(u_ddm3y, r, q, title=f'DDM3Y - ZR - {title}')
print_all(u_m3y, r, q, title=f'M3Y - ZR - {title}')