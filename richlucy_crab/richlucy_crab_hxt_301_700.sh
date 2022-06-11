#!/bin/sh

nskyx=101
nskyy=101
ndetx=80
ndety=80

datafile=../../hitomi/hxt/gimage_hxt1/ah100044010hx1_p1camrec_cl_301_700_x1175y1185_w80.fits
fixed_src_norm_file=skyorg/crab_pulsar_norm.fits
skyfile=none
resp_file=resp/arb_resp.fits
eff_file=resp/arb_eff.fits
nskyx=$nskyx
nskyy=$nskyy
ndetx=$ndetx
ndety=$ndety
outdir=hxt_crab
nem=300
tol_em=1.0e-7
npm=1000
tol_pm=1.0e-7
nnewton=100
tol_newton=1.0e-5

mu_list='1e0 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10' 
for mu_this in $mu_list
do
    mu=$mu_this
    outfile_head=hxt_301_700_${mu}
    ~/work/github/moriiism/srt/richlucy_crab/richlucy_crab \
	$datafile \
	$fixed_src_norm_file \
	$skyfile \
	$resp_file \
	$eff_file \
	$nskyx \
	$nskyy \
	$ndetx \
	$ndety \
	$outdir \
	$outfile_head \
	$nem \
	$tol_em \
	$npm \
	$tol_pm \
	$nnewton \
	$tol_newton \
	$mu
done
