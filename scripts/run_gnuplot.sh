#!/bin/bash

SNPS=$1

cat /dev/null > bar.dat

# From 1 IPU-M2000 (IPU-POD4) to 16 IPU-M2000 (IPU-POD64) accelerators
for POD in 4 8 16 32 64
do
	echo "POD${POD} " >> bar.dat

	for SAMPLES in 4000 8000

	do
		truncate -s -1 bar.dat
		echo -n " " >> bar.dat
		cat output/pod${POD}_${SNPS}snps_${SAMPLES}samples.txt | awk -F ": " '/Tera unique sets per sec./ {print $2}' >> bar.dat
	done
done

gnuplot scripts/gnuplot_drawchart.p
mv output_chart.png output_chart_${SNPS}snps.png
