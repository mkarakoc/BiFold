## BiFold: A Python code for the calculation of double folded (bifold) <br> potentials with density-in/dependent nucleon-nucleon interactions

BiFold calculates the density-dependent (DDM3Yn, BDM3Yn, CDM3Yn) or independent double folded potentials between two colliding spherical nuclei. It is written in a Python package form to give the ability to use the potentials directly in a nuclear reaction/structure code. In addition to using Woods-Saxon/Fermi or Gaussian functions, the code also allows for the definition of nuclear matter densities using pre-calculated densities in a data file. 

### Install using pip from pypi.org
```
pip install bifold
```

### Local install
First download the *'current'* folder or clone the BiFold repository. And, open a terminal in *'current'* folder.
Then type to install:
```
pip install .
```

### Example
```Python
from bifold import *

r = mesh(zero, 10, 0.05)  # fm
q = mesh(zero,  3, 0.05)  # fm^-1

e_lab = 141.7 # MeV
a_proj =  4

rho_p = f_2prm_gaussian(r, 0.4229, (1/0.7024)**.5)
rho_t = f_2prm_fermi(r, 0.169, 3.60, 0.523)

u = u_m3y_reid_zr(e_lab, a_proj, rho_p, rho_t, r, q)

title = "a + 40Ca @ Elab = 141.7 MeV using M3Y-Reid/ZR"
print_all(u, r, q, title=title)
plot_potentials(u, r, part="all")
```


