#!/bin/bash
# i=1
# cat df_choices.csv | while read dfc ; do
# echo "df_choice="$i
# jobname=`echo $dfc | cut -f 1 -d '.'`
r_script=BCC_shared_kriging_run.r
# In the call, $1 should be the path of the folder to run it in (e.g. '/mydir/rubisco_run/'), $2 should be the datestamp of when you ran the data_prep file in yyyymmdd followed by -normed (e.g. 20250815-normed), $3 should be the refseq species in format e.g. AgavicWTLoc.

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${jobname}-$1-$2
#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --mem=100GB
#SBATCH --time=240:00:00
#SBATCH --output=%j_${jobname}_out.log

cd \$SLURM_SUBMIT_DIR
cd log
echo \$SLURM_JOB_ID >> joblist
cd ..
module load R
Rscript $r_script $1 $2 $3
cd log
echo \$SLURM_JOB_ID >> joblist_success
cd ..

EOF

# ((i=i+1))
# done
