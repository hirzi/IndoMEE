#!/bin/bash

# Load modules
source /home/hl636/rds/hpc-work/py37/bin/activate

# Calculate SNB - output dissimilarity matrix
./1.Calculate_SNB.py -f SNB_input_1.txt -r species -m --c1 0 --c2 0 -o SNB_matrixOutput_1_.txt
