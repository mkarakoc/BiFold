from bifold import *
#from bifold.filon import *

r = mesh(zero, 10, 0.05)  # fm
q = mesh(zero,  3, 0.02)  # fm^-1
R = r.copy()

z_proj, a_proj = 2, 4
z_targ, a_targ = 2, 4

rho_p = f_2prm_gaussian(r, 0.4229, (1/0.7024)**.5)
rho_t = rho_p.copy()
vnn = f_2prm_gaussian(r, -33.5606, (1/0.5889)**.5)

# analytical empirical function obtanied by Buck et al.
u1 = f_2prm_gaussian(R, -122.6225, (1/0.22)**.5)

# DFPOT calculation
u2 = f_external(R, r'.\a_a_dfpot_FORT4.txt', data_format=2)

# BiFold calculation
u3 = u_bifold_d(rho_p, rho_t, vnn, r, q)
print_all(u3, r, q)

# compare u1, u2 and u3.
import matplotlib.pyplot as plt
markevery = int(len(R)*5/100)

plt.grid(color='lightgray')
plt.plot(R, -u1.value, label='Buck et al.', linestyle='dotted',
         marker='o', markevery=10, color='black', alpha=0.5, lw=2, ms=9)
plt.plot(R, -u2.value, label='dfpot',linestyle = 'solid',
           marker='s', markevery=15, color='blue', alpha=0.7, lw=2, ms=7)
plt.plot(R, -u3['func_r']['u_R'], label='bifold',linestyle = 'dashed',
           marker='d', markevery=25, color='red', alpha=0.8, lw=2, ms=6)

plt.legend(loc='best')
plt.show()
