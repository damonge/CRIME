#!/bin/bash

source config.sh

fg_type=$1
dirname="${observation_dirname}"
fname_params=${dirname}"/param.ini"

cat > ${fname_params} <<EOF
#Paths and file names
prefix_cosmo= ${cosmo_dirname}/${cosmo_prefix}
prefix_synchrotron= ${synchrotron_dirname}/${synchrotron_prefix}
prefix_galactic_freefree= ${galactic_freefree_dirname}/${galactic_freefree_prefix}
prefix_extragalactic_freefree= ${extragalactic_freefree_dirname}/${extragalactic_freefree_prefix}
prefix_point_sources= ${point_sources_dirname}/${point_sources_prefix}
prefix_custom= none
prefix_out= ${observation_dirname}/${observation_prefix}
fname_nutable= ${cosmo_dirname}/${cosmo_prefix}_nuTable.dat
fname_mask= ${mask_filename}

#Run parameters
seed= ${seed}
nside= ${n_side_joint}
do_polarization= ${do_polarization}

#Instrument parameters
dish_diameter= ${dish_diameter}
n_dishes= ${n_dishes}
t_system= ${t_system}
time_total= ${time_integration}
polarization_leakage= ${polarization_leakage}
EOF
