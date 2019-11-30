#!/usr/bin/env python

import numpy as np

# starting with simple two bits cases
# NOTE: although we are not using 2nd prism slip,
#       need some value as a holder
ratios = [[1, 1, 4, 1, 1],
          [1, 1, 4, 1, 2], [1, 1, 4, 2, 1], [1, 2, 4, 1, 1], [2, 1, 4, 1, 1],
          [1, 1, 4, 2, 2], [1, 2, 4, 1, 2], [2, 1, 4, 1, 2], [1, 2, 4, 2, 1],
          [2, 1, 4, 2, 1], [2, 2, 4, 1, 1],
          [1, 2, 4, 2, 2], [2, 1, 4, 2, 2], [2, 2, 4, 1, 2], [2, 2, 4, 2, 1]]
# reference CRSS (100MPa range)
tau0 = 1e8
taus = 2e8

# prepare common phase section
str_top = "[CP_phenopowerlaw]\n"
str_top += '''
elasticity              hooke
plasticity              phenopowerlaw

(output)                resistance_slip
(output)                shearrate_slip
(output)                resolvedstress_slip
(output)                totalshear
(output)                totalvolfrac
(output)                accumulatedshear_slip

lattice_structure       hex
covera_ratio            1.587
Nslip                   3  3  0  6  12    # per family
Ntwin                   0  0  0  0        # per family

c11                     162.2e9
c12                     91.8e9
c13                     68.8e9
c33                     180.5e9
c44                     46.7e9            # from Webelement

gdot0_slip              0.001             # reference shear rate
n_slip                  50                # might need some calibration
'''
str_bot = '''
a_slip                  2.6454E-01
gdot0_twin              0.001   # no effect
n_twin                  50
tau0_twin               1.06e9  1.06e9  1.06e9  1.06e9
s_pr                    5.88e7
twin_b                  2
twin_c                  25
twin_d                  0.1
twin_e                  0.1
h0_slipslip             4.8703E+08
h0_sliptwin             0.0                # twinning is not active
h0_twinslip             7.38e8
h0_twintwin             4.70e8

# self hardening is set to 1.0 and latent hardening is set to 1.4
interaction_slipslip    1.0 1.4 1.0 1.0 1.4 1.0 1.0 1.0 1.0 1.4 1.0 1.0 1.0 1.0 1.0 1.0 1.4 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.4 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.4 1.0 1.0 1.0 1.0 1.0
interaction_sliptwin    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0                      # no effect
interaction_twinslip    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1                      # no effect
interaction_twintwin    1 1 1 1 1 1 1 1 10 10 10 10 10 10 10 10 10 10 10 10  # no effect
'''

for r in ratios:
    tau0s = np.array(r) * tau0
    tauss = np.array(r) * taus
    outstr = str_top + "tau0_slip\t" + "\t".join(map(str, tau0s)) + "\n"
    outstr += "tausat_slip\t" + "\t".join(map(str, tauss))
    outstr += str_bot
    # configure file name: phase.ratios.config
    fname = "phase." + "-".join(map(str, r)) + ".config"
    with open(fname, 'w') as f:
        f.write(outstr)
