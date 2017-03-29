#!/bin/bash -l
################################################################################
#	SET-UP
################################################################################
echo "SINGLE STRANDED DNA ANALYSIS"
read -r -p "Please write the file name: "  sys
read -r -p "Simulation?: " sim

RED='\033[0;31m'
BLACK='\033[0m'

t0="0"
stride="1"
ncount="2"

# Prepare folders
mkdir -p $(pwd)/${sys}/Analysis/
mkdir -p $(pwd)/${sys}/Data/
mkdir -p $(pwd)/${sys}/Data/${sim}
folder=$(pwd)/${sys}

# Copy raw code files
cp ./Codes/settime.py ${folder}/Simulations/${sim}/
gfortran ./Codes/_index.f08 -o ${folder}/Simulations/${sim}/index
gfortran ./Codes/_prm_calc.f08 -o ${folder}/Simulations/${sim}/prm_calc

# Change to working directory 
cd ${folder}/Simulations/${sim}

# Prepare environment files
echo $sys > system.inp
echo $ncount >> system.inp

# Set topology and trajectory files as variables
top="${folder}/Simulations/${sim}/${sys}_0.tpr"
xtc="${folder}/Simulations/${sim}/${sys}.xtc"

################################################################################
#	MAKE NDX FILES
################################################################################

clear;clear;echo ${RED}"Creating Index Files"${BLACK}
# 1) DNA Index File
echo "keep 1" > ${sys}.inp
echo "q" >> ${sys}.inp
gmx_mpi make_ndx -f ${top} -o ${sys}.ndx < ${sys}.inp
clear;clear;echo "Index: 1/2 COMPLETE"

# 2) Nucleic Base Atomic Site Index Files
./index;rm -f index
clear;clear;echo "Index: 2/2 COMPLETE"

################################################################################
#	TRAJECTORY CORRECTIONS
################################################################################
# The following creates time dependent coordinates in a ".g96" format, which is
# more easily read by FORTRAN.

#clear;clear;echo ${RED}"Correcting trajectory files"${BLACK}
echo "1" > trjconv.inp
# 1) DNA
python settime.py 
gmx_mpi trjcat -f ${sys}_*0.xtc -o ${xtc} -b ${t0} -settime < settime.inp
gmx_mpi trjconv -quiet -f ${xtc} -s ${top} -o ${xtc} -pbc cluster -n ${sys}.ndx
echo "Periodic Boundary Conditions"
gmx_mpi trjconv -quiet -f ${xtc} -s ${top} -o ${xtc} -center -n ${sys}.ndxx
gmx_mpi trjconv -quiet -f ${xtc} -s ${top} -o ${xtc} -fit rot+trans -n ${sys}.ndx
echo "Rotational and Translational Diffussion"

# 2) Nucleic Base Atomic Sites
gmx_mpi trjconv -quiet -f $xtc -s $top -n Atoms.ndx -o Atoms.g96 -dt $stride

#################################################################################
##	CALCULATE TIME DEPENDENT PARAMETER
#################################################################################

echo ${RED}"Calculating structural parameters"${BLACK}
./prm_calc;rm prm_calc
mv Distance ${folder}/Data/${sim}/
mv Roll ${folder}/Data/${sim}/
mv Tilt ${folder}/Data/${sim}/
mv Twist ${folder}/Data/${sim}/

################################################################################
#	CLEAN-UP
################################################################################

#rm -f *.inp *.f08 *.g96 
cd ${folder}/Simulations/${sim}/;rm -f *.ndx *#*#
cd ${folder}

