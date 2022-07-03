#from bifold import *
from bifold.filon import *

r = mesh(zero, 10, 0.05)  # fm
q = mesh(zero,  3, 0.02)  # fm^-1
R = r.copy()
s = r.copy()

# Phys. Rev. C 56 (1997) 954.
# Dao T. Khoa, G. R. Satchler, and W. von Oertzen
# Fig 5. in page 961.

# Physics Letters B 342(1995) 6 - 12
# D.T. Khoa, W. von Oertzen.
# Fig 3. in page 8.
u_fig5_bdm3y1_t_poly = [-0.000136496, 0.00570176, -0.0862306, 0.387249, -0.476532, 4.03722, -0.200097, -117.938]
u_fig5_bdm3y1_t = f_extrapolate(R, 0.01254, 5.04892, u_fig5_bdm3y1_t_poly)
plot_fig5 = pplot(R, u_fig5_bdm3y1_t, label='total - Fig 5.',linestyle = 'None', marker='s',
                  markevery=int(len(R)*2.5/100), color='green', alpha=0.5)

e_lab = 104.0
z_proj, a_proj =  2,  4
z_targ, a_targ = 20, 40
Cs = '1/4'

rho_p = f_2prm_gaussian(r, 0.4229, (1/0.7024)**.5)
rho_t = f_2prm_fermi(r, 0.169, 3.60, 0.523)

rc = 1.405 * (power(a_proj, 1 / 3) + power(a_targ, 1 / 3))
u_coul = u_coul_ucs(R, rc, z_proj, z_targ)

# de Vries
Ri = [0.2,0.6,0.9,1.4,1.9,2.3,2.6,3.1,3.5,4.2,4.9,5.2]
Qi = [0.034724,0.430761,0.203166,0.192986,0.083866,
      0.033007,0.014201,0.000000,0.006860,0.000000,
      0.000438,0.000000]
RP = 1 # fm
rho_p_ch = f_sog(r, Ri, Qi, RP, Ze = 1.9999820)
rho_t_ch = f_3prm_fermi(r, 0.084933, -0.161/3.766**2, 3.766, 0.586)

v_coul = v_coulomb(r)
u_coul_df = u_coul_bifold_d(rho_p_ch, rho_t_ch, v_coul, r, q, R, s)
print_all(u_coul_df, r, q, title='Coulomb potential using double folding model.')
plot_potentials([u_coul, u_coul_df], R)


u_bdm3y1 = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul,
                  r, q, R, s, Cs=eval(Cs), dd_name='bdm3y1', vnn_name='paris')
u_bdm3y1_coul_df = u_xdm3yn_fr(e_lab, a_proj, a_targ, rho_p, rho_t, u_coul_df,
                  r, q, R, s, Cs=eval(Cs), dd_name='bdm3y1', vnn_name='paris')
print_all(u_bdm3y1, r, q, title='exchange part using analytical Coulomb potential')
print_all(u_bdm3y1_coul_df, r, q, title='exchange part using folded Coulomb potential')


suptitle = f'BiFold vs literature for BDM3Y1/Paris, Cs = {Cs}'
title ='Reference: Fig. 5. @ Phys. Rev. C 56 (1997) 954.'
plot_potentials([u_bdm3y1, u_bdm3y1_coul_df], R, part='all', legend_ext=[' - coul_ucs', ' - coul_df'],
                add_plot=plot_fig5, xlimit=(0,5), title=title, suptitle=suptitle)
