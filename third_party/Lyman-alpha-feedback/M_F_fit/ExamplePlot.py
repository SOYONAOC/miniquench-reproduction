""" Some example plots of M_F, as predicted by the fitting function. """

from M_F_fit import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.font_manager as font_manager
import warnings
warnings.filterwarnings("ignore")

###############################################################
#######                                                 #######
#######           P L O T   C O S M E T I C S           #######
#######                                                 #######
###############################################################

# COMMENT THIS OUT IF YOU DON'T HAVE THE REQUIRED FONTS:

#plt.rc('font',**{'sans-serif':['RomanD'],
                 #'size': 27} )
#plt.rc('font',**{'sans-serif':['AVHershey Simplex'],
                #'size': 27} )
plt.rc('font',**{'sans-serif':['AVHershey Complex'],
                'size': 37} )
#plt.rc('font',**{'sans-serif':['Complex'],
                # 'size': 27} )

params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)

plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'GreekC'

# The stuff below need not be removed:

params = {'legend.fontsize': 19 }
plt.rcParams.update(params)

plt.rcParams["legend.markerscale"] = 0.75

plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.major.width'] = 1.7
plt.rcParams['xtick.minor.width'] = 1.7
plt.rcParams['xtick.major.size'] = 20
plt.rcParams['xtick.minor.size'] = 8
plt.rcParams['ytick.major.width'] = 1.7
plt.rcParams['ytick.minor.width'] = 1.7
plt.rcParams['ytick.major.size'] = 20
plt.rcParams['ytick.minor.size'] = 8
plt.rcParams['xtick.major.pad']='15'
plt.rcParams['ytick.major.pad']='15'

plt.rcParams['legend.frameon'] = 'False'

plt.rcParams['axes.linewidth'] = 3.2
plt.rcParams['grid.linewidth'] = 3.2
plt.rcParams['grid.linestyle'] = '--'

###############################################################
#######                                                 #######
#######               P L O T  I T !                    #######
#######                                                 #######
###############################################################

fig, ax = plt.subplots(2, 2, sharey = True)
fig.set_size_inches(15.0, 16.0)

N_H  = np.logspace(19,24)
vmax = np.logspace(-2,3)
T    = 100.0

ax[0,0].set_xscale('log'); ax[0,0].set_yscale('log')
ax[0,1].set_xscale('log'); ax[0,1].set_yscale('log')
ax[1,0].set_xscale('log'); ax[1,0].set_yscale('log')
ax[1,1].set_xscale('log'); ax[1,1].set_yscale('log')

# Effect of dust, point source: 

ax[0,0].plot( N_H, M_F_fit(N_H, T, 0.0, 0.0), linewidth = 3.0,
              color = 'dodgerblue' )
ax[0,0].plot( N_H, M_F_fit(N_H, T, 1e-4, 0.0), linewidth = 3.0,
              color = 'deeppink')
ax[0,0].plot( N_H, M_F_fit(N_H, T, 1e-2, 0.0), linewidth = 3.0,
              color = 'purple' )
ax[0,0].plot( N_H, M_F_fit(N_H, T, 1.0, 0.0), linewidth = 3.0,
              color = 'red' )

ax[0,0].text(3e19, 20, r'D/D$_{\odot}$ = 0', color = 'dodgerblue', fontsize = 28)
ax[0,0].text(3e19, 0.5*20, r'D/D$_{\odot}$ = 10$^{-4}$', color = 'deeppink', fontsize = 28)
ax[0,0].text(3e19, (0.5**2)*20, r'D/D$_{\odot}$ = 10$^{-2}$', color = 'purple', fontsize = 28)
ax[0,0].text(3e19, (0.5**3)*20, r'D/D$_{\odot}$ = 1', color = 'red', fontsize = 28)

ax[0,0].text(4e21, 3, 'Point source'+'\n'+'T = 100 K'+'\n'+'Static cloud', fontsize = 28)
ax[0,0].set_ylabel(r'Force multiplier M$_{F}$', labelpad = 20)
ax[0,0].set_xlabel(r'N$_{HI}$ (cm$^{-2}$)', labelpad = 20)
ax[0,0].set_ylim(1.0, 1e4)

yl = [1.0, 10.0, 100.0, 1e3, 1e4]
ylabels = ['1', '10', r'10$^2$', r'10$^3$', r'10$^4$']
ax[0,0].set_yticks(yl)
ax[0,0].set_yticklabels(ylabels)

# Effect of dust, uniform source: 

ax[0,1].plot( N_H, M_F_fit(N_H, T, 0.0, 1.0), linewidth = 3.0,
              color = 'dodgerblue' )
ax[0,1].plot( N_H, M_F_fit(N_H, T, 1e-4, 1.0), linewidth = 3.0,
              color = 'deeppink' )
ax[0,1].plot( N_H, M_F_fit(N_H, T, 1e-2, 1.0), linewidth = 3.0,
              color = 'purple' )
ax[0,1].plot( N_H, M_F_fit(N_H, T, 1.0, 1.0), linewidth = 3.0,
              color = 'red' )

ax[0,1].text(1e21, 8e2, 'Uniform source'+'\n'+'T = 100 K'+'\n'+'Static cloud', fontsize = 28)
ax[0,1].set_xlabel(r'N$_{HI}$ (cm$^{-2}$)', labelpad = 20)
ax[0,1].set_ylim(1.0, 1e4)

ax[0,1].set_yticks(yl)
ax[0,1].set_yticklabels(ylabels)

# Effect of velocity gradients, point source,
# N_HI = 1e22 cm^-2:

ax[1,0].plot( vmax, M_F_fit(1e22, T, 0.0, 0.0, vmax), linewidth = 3.0,
              color = 'dodgerblue' )
ax[1,0].plot( vmax, M_F_fit(1e22, T, 1e-4, 0.0, vmax), linewidth = 3.0,
              color = 'deeppink')
ax[1,0].plot( vmax, M_F_fit(1e22, T, 1e-2, 0.0, vmax), linewidth = 3.0,
              color = 'purple' )
ax[1,0].plot( vmax, M_F_fit(1e22, T, 1.0, 0.0, vmax), linewidth = 3.0,
              color = 'red' )

ax[1,0].text(1.0, 3.0, 'Point source'+'\n'+'T = 100 K'+'\n'+r'N$_{HI}$ = 10$^{22}$ cm$^{-2}$',
             fontsize = 28)
ax[1,0].set_ylabel(r'Force multiplier M$_{F}$', labelpad = 20)
ax[1,0].set_xlabel(r'$\dot{R}_{cl}$ (km s$^{-1}$)', labelpad = 20)
ax[1,0].set_ylim(1.0, 1e4)

ax[1,0].set_yticks(yl)
ax[1,0].set_yticklabels(ylabels)
xl = [0.01, 0.1, 1.0, 10.0, 100.0, 1e3]
xlabels = ['0.01', '0.1', '1', '10', r'10$^2$', r'10$^{3}$']
ax[1,0].set_xticks(xl)
ax[1,0].set_xticklabels(xlabels)

# Effect of velocity gradients, uniform source,
# N_HI = 1e22 cm^-2:

ax[1,1].plot( vmax, M_F_fit(1e22, T, 0.0, 1.0, vmax), linewidth = 3.0,
              color = 'dodgerblue' )
ax[1,1].plot( vmax, M_F_fit(1e22, T, 1e-4, 1.0, vmax), linewidth = 3.0,
              color = 'deeppink')
ax[1,1].plot( vmax, M_F_fit(1e22, T, 1e-2, 1.0, vmax), linewidth = 3.0,
              color = 'purple' )
ax[1,1].plot( vmax, M_F_fit(1e22, T, 1.0, 1.0, vmax), linewidth = 3.0,
              color = 'red' )

ax[1,1].text(1.0, 400.0, 'Uniform source'+'\n'+'T = 100 K'+'\n'+r'N$_{HI}$ = 10$^{22}$ cm$^{-2}$',
             fontsize = 28)
ax[1,1].set_xlabel(r'$\dot{R}_{cl}$ (km s$^{-1}$)', labelpad = 20)
ax[1,1].set_ylim(1.0, 1e4)

ax[1,1].set_yticks(yl)
ax[1,1].set_yticklabels(ylabels)
ax[1,1].set_xticks(xl)
ax[1,1].set_xticklabels(xlabels)

plt.tight_layout(pad = 0.3)
plt.subplots_adjust(hspace=0.4)
plt.savefig('MF_exampleplot.png', dpi = 300)
plt.show()
