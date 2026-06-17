#!/bin/bash

python3 data_prep.py

datestamp=$(date +"%Y%m%d"-normed)
project="rubisco"
r_script=BCC_shared_kriging_run.r

# module load R
for item in Agavic Aratha Atrund Caredu Cerdem Chealb Chlrei Crimar Cycpan Dacglo Drocap Droven Echele Euggra Flabid Flaflo Flapri Iridou Limlat Limste Limvir Marpol Metgly Nictab Nymalb Palmar Plasup Pteaqu Sarfla Sidcre Spiole Teuhet Trafor Triaes Zeamay;
#for item in Agavic;
    do Rscript $r_script ${project} ${datestamp} ${item}WTLoc; 
done
