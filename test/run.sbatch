#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00

. "/homes/liu3zhen/anaconda3/etc/profile.d/conda.sh"
export PATH="/homes/liu3zhen/anaconda3/bin:$PATH"
conda activate syri2

qry_fas=data/qry.fasta
ref_fas=data/ref.fasta

qry_chr=chr7
ref_chr=chr7

qry_name=Q
ref_name=R

perl ~/scripts2/chrcomp/chrcomp \
	--qry $qry_fas \
	--qchr $qry_chr \
	--qryname $qry_name \
	--ref $ref_fas \
	--rchr $ref_chr \
	--refname $ref_name \
	--match 10000 \
	--identity 95 \
	--minsyn 10000 \
	--minothers 50000

