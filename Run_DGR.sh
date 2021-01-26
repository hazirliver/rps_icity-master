#!/bin/bash

mkdir RT
mkdir TR
mkdir VR

python3 main.py -s seeds_rt.tsv -f _RT
mv *_RT* ./RT

python3 main.py -s seeds_tr.tsv -f _TR
mv *_TR* ./TR

python3 main.py -s seeds_vr.tsv -f _VR
mv *_VR* ./VR

mkdir Venn
Rscript Venn.R