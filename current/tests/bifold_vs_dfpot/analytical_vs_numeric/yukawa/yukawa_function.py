from bifold import *

def u_analytically_folded(r, rho0, a, v0, c):
    u_numerator = 4*pi*(48*a**5*exp(a*r)*pi + exp(c*r)*pi*(-48*a**5 -
                        3*(a - c)*(a + c)*(11*a**4 - 4*a**2*c**2 + c**4)*r -
                        3*a*(a**2 - c**2)**2*(3*a**2 - c**2)*r**2 -
                        a**2*(a**2 - c**2)**3*r**3))*v0*rho0**2

    u_denominator = 3*a**3*c*(a**2 - c**2)**4*exp((a + c)*r)*r
    u = u_numerator/u_denominator
    u[0] = f_0(r, u)
    return u

rho0, a = 4.03242, 4.66232
v1, b1, v2, b2 = +7999.00, 4.0, -2134.25, 2.5

r = mesh(zero,  5, 0.02)  # fm
q1 = mesh(zero,10, 0.02)  # fm^-1
q2 = mesh(zero, 3, 0.02)  # fm^-1

rho_p = f_exp_decay(r, rho0, a)
rho_t = rho_p.copy()
vnn = v_m3y_reid_d(r)
u_bf1 = u_bifold_d(rho_p, rho_t, vnn, r, q1)
u_bf2 = u_bifold_d(rho_p, rho_t, vnn, r, q2)

u_f1 = u_analytically_folded(r, rho0, a, v1, b1)
u_f2 = u_analytically_folded(r, rho0, a, v2, b2)

print_all(u_bf1, r, q1)
print_all(u_bf2, r, q2)

analytic = pplot(r, u_f1 + u_f2, label='analytic',linestyle = 'None',
                 marker='o', markevery=10, color='red', alpha=0.5, lw=2, ms=6)
plot_potentials([u_bf1, u_bf2], r, add_plot=analytic, colors=['black', 'blue','orange'],
                legends=['bifold - $q_{max}$ = 10 fm$^{-1}$', 'bifold - $q_{max}$ = 3 fm$^{-1}$'])