#!/bin/bash

# Loads specific versions of modules 
#module load graphcore/vipu/1.17.1
#module load graphcore/sdk/2.6.0

mkdir bin
mkdir datasets
mkdir output

# Compiles binaries targeting different amounts of IPU processors
make

# Generates datasets with required dimensions
sh scripts/generate_datasets.sh $1

export IPUOF_CONFIG_PATH=/cm/shared/apps/graphcore/vipu/etc/ipuof.conf.d/p64_cl_a01_a16.conf
sh scripts/run_pods.sh $1

# Generates plot for visualization
sh scripts/run_gnuplot.sh $1
