from bifold import *
#from bifold.filon import *

r = mesh(zero, 10, 0.05)  # fm
q = mesh(zero,  3, 0.02)  # fm^-1
R = r.copy()
s = r.copy()

# Phys. Rev. C 56 (1997) 954.
# Dao T. Khoa, G. R. Satchler, and W. von Oertzen
# Fig 5. in page 961.

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