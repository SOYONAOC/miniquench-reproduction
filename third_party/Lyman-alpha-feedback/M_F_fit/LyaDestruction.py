""" This module contains everything related to Lya continuum absorption
    and the Lya destruction probability. """

import numpy as np
from ConstantsParameters import *

###############################################################
#######                                                 #######
#######      C O N T I N U U M  A B S O R P T I O N     #######
#######                                                 #######
###############################################################


def tau_cabs(D_Dsun, N_H, N_H2, sigma_dHsun):
    """ The cloud continuum absorption optical depth.
        Arguments:

        D_Dsun:      Dust-to-gas ratio in Solar units
        N_H:         Average total hydrogen (HI+HII+H2) column density (cm^-2)
        N_H2:        Average H2 column density (cm^-2)
        sigma_dHsun: Dust absorption cross-section at Lya at Solar
                     dust abundance, per H nucelon (cm^2/H) """

    # Dust absorption optical depth:

    tau_dabs  = 100.0*(sigma_dHsun*D_Dsun/1e-21)*(N_H/1e23)

    # Raman scattering optical depth:

    tau_H2abs = 2.1e-4*(N_H2/1e21)

    # Total continuum absorption optical depth:

    tau_cabs  = tau_dabs + tau_H2abs

    return tau_cabs

###############################################################
#######                                                 #######
#######  D E S T R U C T I O N  P R O B A B I L I T Y   #######
#######                                                 #######
###############################################################

def p_d(nHI, nHeI, nH2, T, z, M = 0.0, b_s = 0.4):
    """ The Lya destruction probability. Arguments:

        nHI:   The average H I number density (cm^-3)
        nHeI:  The average He I number density (cm^-3)
        nH2:   The average H2 number density (cm^-3)
        T:     Gas temperature
        z:     Redshift
        M:     3D turbulent Mach number (controls turbulent fluctuations)
        b_s:   Turbulent driving parameter (1/3 < b_s < 1). (default = 0.4, i.e. 'natural' mix) """

    # CMB stimulated transition rates:

    Gamma_CMB2p2s = 3.1e-6*(1.0+z)
    Gamma_CMB2s2p = 3.0*Gamma_CMB2p2s

    # 2s -> 2p rates for He I, H2, and H I:

    T_100  = T/100.0         # Temperature normalization
    Turb   = 1 + (b_s*M)**2  # Clumping factor for turbulent medium

    C_HeI  = 2.6e-9*(T_100**0.315)*nHeI*Turb
    
    C_H2   = 1.7e-9*(T_100**0.3)*nH2*Turb
    
    C_HI   = ( 4.1e-11*np.exp(-(1.6/T)**(3/4) - 0.15*T)
           + (1.6e-10 + 1.57e-12*(T**1.1))/(
            (8.9*(T**(-0.55)) + 0.078*(T**0.35))**2.0))*nHI*Turb

    # Total 2s <-> 2p rates:

    C_2s2p = C_HeI + C_H2 + C_HI
    C_2p2s = C_2s2p/3.0

    # Destruction probability from 2p -> 2s:

    R_2p2s = (2/3)*A_2p2s + C_2p2s + Gamma_CMB2p2s
    R_2s2p = A_2s2p + C_2s2p + Gamma_CMB2s2p

    p_dHI  = R_2p2s/((A_Lya/A_2ph)*(A_2ph + R_2s2p) + R_2p2s)

    # Destruction probability from H2 line absorption:

    p_dH2  = np.minimum((33.0*np.exp(-14540.0/T) + 13.0*np.exp(-15564.0/T))/(
             46.0 + 1e3*np.exp(-(1e3/T)**1.5)), 1.3e-3)*(nH2/nHI)

    # Total destruction probability:

    p_d    = p_dHI + p_dH2

    return p_d
