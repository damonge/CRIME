#!/bin/bash

source config.sh

fname_params="${cosmo_dirname}/param.ini"
prefix_out="${cosmo_dirname}/${cosmo_prefix}"

cat > ${fname_params} <<EOF
#File names
#Output prefix. Output will be in prefix_###.fits
prefix_out= ${prefix_out}
#Path to power spectrum at z=0. Power spectrum file must
#be in CAMB format: k (h/Mpc), P(k) (Mpc/h)^3.
pk_filename= ${pk_filename}

#Cosmological parameters
#Non-relativistic matter
omega_M= ${omega_M}
#Dark energy
omega_L= ${omega_L}
#Baryons
omega_B= ${omega_B}
#Hubble parameter (in units of 100 km/s/Mpc)
h= ${hhubble}
#Dark energy equation of state
w= ${weos}
#Primordial scalar spectral index
ns= ${nscal}
#Power spectrum normalization
sigma_8= ${sigma_8}

###Extra Gaussian smoothing scale [Mpc/h] (set to a
#negative value if you don't want any smoothing)
r_smooth= ${r_smooth}

### Frequency range for output
#File containing frequency channels
frequencies_filename= ${frequencies_filename}

#Angular resolution for output
n_side= ${n_side}

#Grid parameters
#Will use a Cartesian grid with n_grid^3 cells
n_grid= ${n_grid}

#RNG seed
seed= ${seed}

#Point sources
do_psources= 0
EOF
