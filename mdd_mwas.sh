#! /usr/bin/bash
#$ -l h_rt=0:30:00
#$ -l h_vmem=192G
#$ -cwd
#$ -o logs
#$ -e logs
#$ -b y

script_dir=$(dirname $0)

. /etc/profile.d/modules.sh

module load igmm/apps/R/3.4.1

Rscript $script_dir/mdd_mwas.R $*
