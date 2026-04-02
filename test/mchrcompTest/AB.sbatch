#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00

. "/homes/liu3zhen/anaconda3/etc/profile.d/conda.sh"
export PATH="/homes/liu3zhen/anaconda/bin:$PATH"
conda activate chrcomp26

perl ../../mchrcomp \
	--qry ../data/ATCC64557.fasta \
	--ref ../data/B71.fasta \
	--comp ../data/ABcompTable \
	--prefix outAB
	--match 10000 \
	--identity 95 \
	--minsyn 10000 \
	--minothers 5000

