#!/bin/bash

datestamp="20250825-normed"
rundir=""

for item in Agavic Aratha Atrund Caredu Cerdem Chealb Chlrei Crimar Cycpan Dacglo Drocap Droven Echele Euggra Flabid Flaflo Flapri Iridou Limlat Limste Limvir Marpol Metgly Nictab Nymalb Palmar Plasup Pteaqu Sarfla Sidcre Spiole Teuhet Trafor Triaes Zeamay;
#for item in Agavic;
    do sbatch ./BCC_run_script.sh ${rundir} ${datestamp} ${item}WTLoc; 
done
