#!/bin/bash

source config.sh

echo "Starting... on ${name} ${seed}" $(date)

#Cosmological signal
echo "Running cosmological signal" $(date)
mkdir -p ${cosmo_dirname}
./write_GetHI_params.sh
${GETHI} ${cosmo_dirname}/param.ini
echo " "

#Galactic synchrotron
echo "Running galactic synchrotron" $(date)
mkdir -p ${synchrotron_dirname}
./write_ForGet_params.sh galactic_synchrotron
${FORGET} ${synchrotron_dirname}/param.ini
echo " "

#Galactic free-free
echo "Running galactic free-free" $(date)
mkdir -p ${galactic_freefree_dirname}
./write_ForGet_params.sh galactic_freefree
${FORGET} ${galactic_freefree_dirname}/param.ini
echo " "

#Extragalactic free-free
echo "Running extragalactic free-free" $(date)
mkdir -p ${extragalactic_freefree_dirname}
./write_ForGet_params.sh extragalactic_freefree
${FORGET} ${extragalactic_freefree_dirname}/param.ini
echo " "

#Point sources
echo "Running point sources" $(date)
mkdir -p ${point_sources_dirname}
./write_ForGet_params.sh point_sources
${FORGET} ${point_sources_dirname}/param.ini
echo " "

#Get observation
echo "Compiling observation" $(date)
mkdir -p ${observation_dirname}
./write_JoinT_params.sh
${JOINT} ${observation_dirname}/param.ini
echo " "

echo "Done!" $(date)
