# to use simpson integration by default
from bifold import *

# to use simpson integration
# from bifold.simpson import *

# to use filon integration
# from bifold.filon import *

# meshes for integrations.
# r is for the densities of nuclei.
# q is for all momentum space representations of fourier transforms.
# R defines separation radius between two nuclei.
# s defines separation between two nucleons.
# usually better to use r for R and s if unless not necessary otherwise.
r = mesh(zero, 10, 0.05)  # fm
q = mesh(zero,  3, 0.02)  # fm^-1
# this is to use same mesh [r] for R and s.
R = r.copy()
s = r.copy()

# e_lab is nuclear bombarding energy in the laboratory reference frame.
# z_ and a_ defines atomic number and atomic mass number, respectively.
e_lab = 141.7
z_proj, a_proj =  2,  4
z_targ, a_targ = 20, 40

# nuclear matter densities for projectile and target nuclei
rho_p = f_2prm_gaussian(r, 0.4229, (1/0.7024)**.5)
rho_t = f_2prm_fermi(r, 0.169, 3.60, 0.523)

# double folded calculation
title_u1 = "alpha + 40Ca @ Elab = 141.7 MeV using M3Y-Reid/ZR"
u1 =  u_m3y_reid_zr(e_lab, a_proj, rho_p, rho_t, r, q, R, s)

# double folded calculation
title_u2 = "alpha + 40Ca @ Elab = 141.7 MeV using M3Y-Paris/ZR"
u2 =  u_m3y_paris_zr(e_lab, a_proj, rho_p, rho_t, r, q, R, s)

# print calculations of BiFold
print_all(u1, R, q, title=title_u1)
print('\n\n\n')
print_all(u2, R, q, title=title_u2)

# calculations obtained using DFPOT for M3Y - Reid [direct part] and zero range [exchange part]
u1_dfpot = f_external(R, r'.\dfpot_calc\reid_m3y_zr_FORT4.txt', data_format=2)
u1_dfpot_d = f_external(R, r'.\dfpot_calc\reid_m3y_zr_direct_FORT4.txt', data_format=2)
u1_dfpot_e = f_external(R, r'.\dfpot_calc\reid_m3y_zr_exchange_FORT4.txt', data_format=2)

# calculations obtained using DFPOT for M3Y - Paris [direct part] and zero range [exchange part]
u2_dfpot = f_external(R, r'.\dfpot_calc\paris_m3y_zr_FORT4.txt', data_format=2)
u2_dfpot_d = f_external(R, r'.\dfpot_calc\paris_m3y_zr_direct_FORT4.txt', data_format=2)
u2_dfpot_e = f_external(R, r'.\dfpot_calc\paris_m3y_zr_exchange_FORT4.txt', data_format=2)

# following line to compare BiFold and DFPOT calculations using matplotlib.
import matplotlib.pyplot as plt
markevery = int(len(R)*5/100)

f, ax = plt.subplots(3, 1, sharex='all', figsize=(6,9))
ax[0].grid(color='lightgray', zorder=-999)
y_loc = (min(u2['func_r']['direct']['u_R']) + max(u2['func_r']['direct']['u_R']))/2
ax[0].text(5,y_loc, 'direct')

ax[0].plot(R, u1['func_r']['direct']['u_R'], linestyle='solid', color='black')
ax[0].plot(R, u1_dfpot_d.value, label='dfpot-reid',linestyle = 'None',
           marker='s', markevery=markevery, color='orange', alpha=0.5)
ax[0].plot(R, u2['func_r']['direct']['u_R'], linestyle='dashed', color='red')
ax[0].plot(R, u2_dfpot_d.value, label='dfpot-paris',linestyle = 'None',
           marker='d', markevery=markevery, color='blue', alpha=0.5)

ax[1].grid(color='lightgray', zorder=-999)
y_loc = (min(u2['func_r']['exchange']['u_R']) + max(u2['func_r']['exchange']['u_R']))/2
ax[1].text(5,y_loc, 'exchange')
ax[1].plot(R, u1['func_r']['exchange']['u_R'], linestyle='solid', color='black', label='bifold-reid')
ax[1].plot(R, u1_dfpot_e.value, label='dfpot-reid',linestyle = 'None',
           marker='s', markevery=markevery, color='orange', alpha=0.5)
ax[1].plot(R, u2['func_r']['exchange']['u_R'], linestyle='dashed', color='red', label='bifold-paris')
ax[1].plot(R, u2_dfpot_e.value, label='dfpot-paris',linestyle = 'None',
           marker='d', markevery=markevery, color='blue', alpha=0.5)
ax[1].legend(loc='center right')

ax[2].grid(color='lightgray', zorder=-999)
y_loc = (min(u1['func_r']['total']['u_R']) + max(u1['func_r']['total']['u_R']))/2
ax[2].text(5, y_loc, 'total')
ax[2].plot(R, u1['func_r']['total']['u_R'], linestyle='solid', color='black')
ax[2].plot(R, u1_dfpot.value, label='dfpot-reid',linestyle = 'None',
           marker='s', markevery=markevery, color='orange', alpha=0.5)
ax[2].plot(R, u2['func_r']['total']['u_R'], linestyle='dashed', color='red')
ax[2].plot(R, u2_dfpot.value, label='dfpot-paris',linestyle = 'None',
           marker='d', markevery=markevery, color='blue', alpha=0.5)

f.subplots_adjust(hspace=0.0)
ax[2].set_xlabel('R [fm]', fontsize=14)
ax[1].set_ylabel('U [MeV]', fontsize=14)
plt.xlim(R[0],R[-1])
plt.show()
