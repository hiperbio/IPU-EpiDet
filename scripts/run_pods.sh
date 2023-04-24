#!/bin/bash

SNPS=$1
# From 1 IPU-M2000 (IPU-POD4) to 16 IPU-M2000 (IPU-POD64) compute engines
for POD in 4 8 16 32 64
do
	for SAMPLES in 4000 8000
	do
		./bin/ipu-epidet_pod${POD} datasets/${SNPS}snps_${SAMPLES}samples.csv > output/pod${POD}_${SNPS}snps_${SAMPLES}samples.txt
	done
done

