#$ -V
#$ -S /bin/bash
#$ -e ~/log
#$ -o ~/log
#$ -cwd
#$ -j y
#$ -l mem_free=20G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=20G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=23:59:59
#$ -e ~/log
#$ -o ~/log


cd ~/src/scripts
module load CBI
module load r
#Rscript graph_test.R --cores=="$cores"
Rscript fig2_graphtest.R --cores=="$cores"
