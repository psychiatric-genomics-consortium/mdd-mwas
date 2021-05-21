#$ -N stradl_ewas_chr
#$ -l h_rt=2:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 6
#$ -cwd
#$ -e logs
#$ -o logs

. /etc/profile.d/modules.sh
module load igmm/apps/R/3.4.1

split_chr_probes_r=$1
mvalues_file=$2
out_dir=$3

Rscript $split_chr_probes_r $mvalues_file $out_dir
