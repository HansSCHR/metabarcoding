#!/bin/bash
#!/usr/bin/python3.6
#$ -q bigmem.q
#$ -N merge
#$ -M hans.schrieke@cirad.fr
#$ -pe parallel_smp 15
#$ -l mem_free=6G
#$ -V
#$ -cwd

module load system/conda/5.1.0
#source ~/virtual-envs/illumina-utils-v2.0.2/bin/activate
source activate illumina_utils
python -V

#cd /media/schriekehans/F60E8EF20E8EAAE7/stage/data/runs_fastq_copie_run3_modif_TEST
cd /homedir/schrieke/stage/oligotyping/merge/data_raw_mod/runs_fastq_copie_run3_modif_TEST #folder where are my sequences

  #echo $ini
iu-gen-configs sample_finale.txt -o output

#cd /media/schriekehans/F60E8EF20E8EAAE7/stage/data/runs_fastq_copie_run3_modif_TEST/output
cd /homedir/schrieke/stage/oligotyping/merge/data_raw_mod/runs_fastq_copie_run3_modif_TEST/output #forlder where are my ini files

for ini in *.ini
do
  #echo $ini
  iu-merge-pairs --min-overlap-size 30 --max-num-mismatches 0 --enforce-Q30-check $ini
  echo $ini "done" > state.txt
done
