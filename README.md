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

### Outputs

```bash
═══════════════════════════════╣a + 40Ca @ Elab = 141.7 MeV using M3Y-Reid/ZR╠═══════════════════════════════

density/interaction               L          norm        renorm      vol2        vol4            msr
──────────────────────────────────────────────────────────────────────────────────────────────────────────────
───────total     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
u_R      : u_m3y_reid_zr          0          None        1.000   -59536.370  -982867.456       16.509 
───────direct    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
u_R      : u_direct               0          None        1.000   -23281.699  -486747.684       20.907 
rho_p    : f_2prm_gaussian        0          None        1.000        4.000        8.543        2.136 
rho_t    : f_2prm_fermi           0          None        1.000       39.908      461.090       11.554 
vnn      : f_yukawa               0          None        1.000     1570.558      588.970        0.375 
vnn      : f_yukawa               0          None        1.000    -1716.459    -1647.805        0.960 
───────exchange  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
u_R      : u_exchange_zr          0          None        1.000   -36254.670  -496119.772       13.684 
rho_p    : f_2prm_gaussian        0          None        1.000        4.000        8.543        2.136 
rho_t    : f_2prm_fermi           0          None        1.000       39.908      461.090       11.554 
vnn      : f_dirac_delta          0          None        1.000     -227.114        0.000        0.000 

  R            u_R       	   R            u_R       	   R            u_R       	   R            u_R       	
───────   ───────────────	 ───────   ───────────────	 ───────   ───────────────	 ───────   ───────────────	
  0.000    -224.563778558	   2.500    -157.142066726	   5.000     -33.591208578	   7.500      -1.294809340	
  0.050    -224.538322266	   2.550    -154.527402817	   5.050     -32.017487683	   7.550      -1.195282003	
  0.100    -224.461928239	   2.600    -151.883273803	   5.100     -30.494524998	   7.600      -1.102914652	
  0.150    -224.334521330	   2.650    -149.211853484	   5.150     -29.022094646	   7.650      -1.017241399	
  ...
  ...
  ...
  2.300    -167.264509196	   4.800     -40.395194787	   7.300      -1.774804137	   9.800      -0.024215666	
  2.350    -164.788335884	   4.850     -38.617900517	   7.350      -1.641422081	   9.850      -0.022007393	
  2.400    -162.274567943	   4.900     -36.891404224	   7.400      -1.517343866	   9.900      -0.019983036	
  2.450    -159.725143549	   4.950     -35.215829117	   7.450      -1.401989622	   9.950      -0.018129460	
```



