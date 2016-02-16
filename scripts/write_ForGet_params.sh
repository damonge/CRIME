#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: write_ForGet_params.sh foreground_type"
    exit 1
fi

source config.sh

fg_type=$1
dirname="null"
fname_params="null"
prefix_out="null"

if [ ${fg_type} == "galactic_synchrotron" ]
then
    echo "Doing galactic synchrotron"
    dirname="${synchrotron_dirname}"
    prefix_out=${dirname}"/${synchrotron_prefix}"
elif [ ${fg_type} == "galactic_freefree" ]
then
    dirname="${galactic_freefree_dirname}"
    prefix_out=${dirname}"/${galactic_freefree_prefix}"
elif [ ${fg_type} == "extragalactic_freefree" ]
then
    dirname="${extragalactic_freefree_dirname}"
    prefix_out=${dirname}"/${extragalactic_freefree_prefix}"
elif [ ${fg_type} == "point_sources" ]
then
    dirname="${point_sources_dirname}"
    prefix_out=${dirname}"/${point_sources_prefix}"
else
    echo "Unknown foreground ${fg_type}."
    echo "Possible values: galactic_synchrotron, galactic_freefree, extragalactic_freefree and point_sources"
    exit 1
fi
fname_params=${dirname}"/param.ini"

cat > ${fname_params} <<EOF
#File names
fname_nutable= ${cosmo_dirname}/${cosmo_prefix}_nuTable.dat
fname_haslam= ${haslam_filename}
fname_specin= ${specin_filename}
prefix_out= ${prefix_out}

#Run parameters
lmin= ${lmin}
lmax= ${lmax}
nside= ${n_side}
seed= ${seed}

#Cl model
cl_model= ${fg_type}
amp= 57.0
beta= 1.1
alpha= 2.07
xi= 1.0
lref= 1000
nu_ref= 130.0

#Polarization
do_polarization= ${do_polarization}
xi_polarization= ${xi_polarization}
beta_polarization= ${beta_polarization}
fname_faraday= ${faraday_filename}
EOF
