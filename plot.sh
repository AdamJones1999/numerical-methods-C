#!/bin/bash
Nvars=$2
Nt=$3
python3 shmODE-error-plot.py data/test$1.data $Nvars $Nt
