#!/bin/bash
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

