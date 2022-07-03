from bifold.filon import *

r = mesh(zero, 10, 0.05)  # fm
q = mesh(zero,  3, 0.02)  # fm^-1
R = r.copy()
s = r.copy()

# Phys. Rev. C 49 (1994) 1652.
u_fig1_ddm3y1_t_poly = {
    160: [-0.00127708, 0.0284608, -0.181577, 0.121312, -0.564074, 20.4557, -1.279950, -340.378],
    480: [-0.00129518, 0.0279258, -0.177686, 0.191960, -0.891316, 18.5787, -1.528820, -290.491],
    960: [-0.00151486, 0.0339716, -0.249294, 0.615374, -1.564210, 14.8922, -1.149480, -225.148],
   1440: [-0.00121332, 0.0267145, -0.186182, 0.347161, -0.446128, 9.73655, -0.436533, -172.022]
}

u_fig1_ddm3y1_t = {
    160: f_extrapolate(R,  0.00580, 7.9696, u_fig1_ddm3y1_t_poly[160]),
    480: f_extrapolate(R, -0.00253, 7.8282, u_fig1_ddm3y1_t_poly[480]),
    960: f_extrapolate(R, -0.00646, 7.8494, u_fig1_ddm3y1_t_poly[960]),
   1440: f_extrapolate(R, -0.00225, 7.9272, u_fig1_ddm3y1_t_poly[1440])
}

z_proj, a_proj =  8, 16
z_targ, a_targ =  8, 16

rho_p = f_2prm_fermi(r, 0.181, 2.525, 0.450)
rho_t = f_2prm_fermi(r, 0.181, 2.525, 0.450)

rc = 1.2 * (power(a_proj, 1 / 3) + power(a_targ, 1 / 3))

u_coul = u_coul_ucs(R, rc, z_proj, z_targ)

e_labs = [160.0, 480, 960, 1440]
u_ddm3y1 = dict()
for e_lab in e_labs:
    u_ddm3y1[e_lab] = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                        r, q, R, s, Cs=1/4, dd_name='ddm3y1', vnn_name='reid')
    print_all(u_ddm3y1[e_lab], r, q, title=f'O^16 + O^16 @ Elab = {e_lab} MeV')


import matplotlib.pyplot as plt
markevery = int(len(R)*2.5/100)
linestyles = ['solid', 'dotted', 'dashed', 'dashdot']
markers = ['s', 'D', '^', 'v']

plt.figure(figsize=(9,6))
plt.suptitle('BiFold vs literature for DDM3Y1/Reid')
plt.title('Reference: Fig. 5. @ Phys. Rev. C 49 (1994) 1652.')
plt.grid(color='lightgray', linestyle='dashed', zorder=-999)

for i, e_lab in enumerate(e_labs):
    plt.plot(R, u_ddm3y1[e_lab]['func_r']['total']['u_R'], label=f'bifold {e_lab} MeV', linestyle=linestyles[i])
    plt.plot(R, u_fig1_ddm3y1_t[e_lab], label=f'total - Fig 5. {e_lab} MeV',linestyle = 'None',
             marker=markers[i], markevery=markevery, alpha=0.5)

plt.legend(loc='lower right')
plt.tick_params(axis='both', which='major', labelsize=12)
plt.xlabel('R [fm]', fontsize=16)
plt.ylabel('U(R) [MeV]', fontsize=16)
plt.xlim(0, 8)
plt.show()
