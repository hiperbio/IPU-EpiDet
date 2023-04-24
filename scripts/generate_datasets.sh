#!/bin/bash

for SAMPLES in 4000 8000
do
	python scripts/dataset_generator.py ${1} ${SAMPLES} > datasets/${1}snps_${SAMPLES}samples.csv
done
