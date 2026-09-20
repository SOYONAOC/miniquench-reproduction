""" Some useful physical constants, parameters, and conversion factors. """

import numpy as np

# Fundamental constants:

c                   =  2.9979245e10       # Speed of ligth (cm/s)
kB                  =  1.3806488e-16      # Boltzmann constant (erg/K)
mH                  =  1.6735327e-24      # Mass of hydrogen atom (g)

# Hydrogen, helium, and dust abundance:

Y                   =  0.247              # Helium mass fraction
X                   =  1.0-Y              # Hydrogen mass fraction
Dsun                =  1/162              # Solar dust-to-gas ratio

# Dust absorption cross-section at Lya at Solar dust abundance,
# per H nucelon (cm^2/H), for different dust models:

sigma_dHsun_HD23    =  6.873e4*Dsun*mH/X  # Hensley & Draine (2023), Astrodust + PAH
sigma_dHsun_WD_Rv3  =  6.854e4*Dsun*mH/X  # Weingartner & Draine (2001), Milky Way, R_V = 3.1
sigma_dHsun_WD_Rv4  =  4.425e4*Dsun*mH/X  # Weingartner & Draine (2001), Milky Way, R_V = 4.0
sigma_dHsun_WD_Rv5  =  3.182e4*Dsun*mH/X  # Weingartner & Draine (2001), Milky Way, R_V = 5.5
sigma_dHsun_WD_LMC  =  7.251e4*Dsun*mH/X  # Weingartner & Draine (2001), LMC average
sigma_dHsun_WD_SMC  =  7.181e4*Dsun*mH/X  # Weingartner & Draine (2001), SMC bar

# Einstein-A coefficients:

A_Lya               =  6.265e8            # Lya Einstein-A coefficient (s^-1)
A_2ph               =  8.2206             # 2s->1s Einstein-A coefficient(s^-1)
A_2p2s              =  8.78e-7            # 2p->2s Einstein-A coefficient(s^-1)
A_2s2p              =  1.597e-9           # 2s->2p Einstein-A coefficient(s^-1)


# Unit conversion factors:

kms                 =  1e5                # 1 km/s in cm/s
pc                  =  3.086e18           # 1 pc in cm
AU                  =  1.496e13           # 1 AU in cm
Msun                =  1.989e33           # Solar mass in g
