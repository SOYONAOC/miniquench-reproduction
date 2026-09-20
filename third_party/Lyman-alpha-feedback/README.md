## Repository for project on Lyman-α radiation pressure feedback
<div style="text-align: center;">
    <img src="https://github.com/user-attachments/assets/fa23c18b-6641-4bc1-beed-804950a428bb" alt="Schematic Overview Lya" width="550"/>
</div>

### Purpose of repository

This repository contains files and data related to the following paper on Lyman-α (Lyα) radiation pressure feedback:

[Nebrin et al. (2025), 'Lyman-α feedback prevails at Cosmic Dawn: Implications for the first galaxies, stars, and star clusters', MNRAS, 537, 2, 1646-1687](https://ui.adsabs.harvard.edu/abs/2025MNRAS.537.1646N/abstract).

Radiation pressure from trapped Lyα photons is one of the most important processes that regulate how efficiently the Universe form stars. 

#### Authors & contributors: Olof Nebrin, Aaron Smith, Kevin Lorinc, Johan Hörnquist, Åsa Larsson, Garrelt Mellema, and Sambit K. Giri 

### Description of folders and files

#### M_F_fit:

This folder contains a Python implementation of the fit to the Lyα force multiplier M<sub>F</sub>, as given the paper, as well as related quantities 
(e.g. the destruction probability). An example plot is included too. If you use these fits or results in your work, cite [Nebrin et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025MNRAS.537.1646N/abstract).

#### MCRT_forcemultiplier:

This folder contains the Monte Carlo radiative transfer (MCRT) data for the Lyα force multiplier, as plotted in Fig. 5 in [Nebrin et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025MNRAS.537.1646N/abstract). Each data file is named in a self-explanatory manner with respect to Fig. 5. In each file, the first column is a<sub>v</sub>τ<sub>cl</sub>, and the second column the force multiplier M<sub>F</sub>. If you use this data, cite [Nebrin et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025MNRAS.537.1646N/abstract). 

#### H_H_cross_section:

The file in this folder contains the data for the cross-section for H(2s) + H(1s) -> H(2p) + H(1s). The first column is collision energy (in eV), and the second column the cross-section (in cm<sup>2</sup>). If you use this data, cite [Nebrin et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025MNRAS.537.1646N/abstract), [Hörnquist et al. (2022)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.106.062821), and [Hörnquist et al. (2023)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.108.052811).

### Contact

If you have any questions regarding the paper or the files in this repository, don't hesitate to reach out at: olof.nebrin@astro.su.se


