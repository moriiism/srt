#!/usr/bin/perl

#
# emgist_case2_sim_cv.pl
#
# 2018.08.21  M.Morii
#

use strict;
use warnings;
use POSIX;
no warnings qw(once);

#stop buffering
$| = 1;


$ENV{'LANG'} = "C";

my $mipllib = sprintf("%s/script/mipllib/mipllib.pl",
		      $ENV{'MITOOL'});
require($mipllib);
$mipllib::sevar{'progname'} = "emgist_case2_sim_cv.pl";

# main
{
    printf("#--- %s: START ---#\n", $mipllib::sevar{'progname'});

    my $cmd;
    # -- argument
    my $NARG = 5;
    my($infile, $nevt, $nfold, $npartition, $respdir);
    
    # infile = ../data/20180313/crab/crab_gmodel_his_sel1_1000000.img
    # nevt = 1000000
    # nfold = 10
    # npartition = 10

    printf("# of arg = %d\n", $#ARGV + 1);
    if( $#ARGV == $NARG - 1){
	my $iarg = 0;
	$infile       = $ARGV[$iarg]; $iarg ++;
	$nevt         = $ARGV[$iarg]; $iarg ++;
	$nfold        = $ARGV[$iarg]; $iarg ++;
	$npartition   = $ARGV[$iarg]; $iarg ++;
	$respdir      = $ARGV[$iarg]; $iarg ++;
    }else{
        printf("# of arg must be %d.\n", $NARG);
	&Usage();
    }
    my $prompt_arg = "arg";
    printf("%s: %s: infile        = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $infile);
    printf("%s: %s: nevt          = %d\n", $mipllib::sevar{'progname'}, $prompt_arg, $nevt);
    printf("%s: %s: nfold         = %d\n", $mipllib::sevar{'progname'}, $prompt_arg, $nfold);
    printf("%s: %s: npartition    = %d\n", $mipllib::sevar{'progname'}, $prompt_arg, $npartition);
    printf("%s: %s: respdir       = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $respdir);
    # == argument

    # simulation
    my $rand_seed_sky = 1;
    my $rand_seed_det = 1;
    my $rand_seed_partition = 1;
    my $outdir_simobs = "simobs";
    my $outfile_head_simobs = "simobs";
    $cmd = sprintf("/home/morii/work/github/moriiism/srt/simobs/simobs  " .
		   "%s  %s  " .
		   "%d  %d  %d  %d  " .
		   "%d  %d  " .
		   "%s  %s  ", 
		   $respdir, $infile,
		   $nevt, $rand_seed_sky, $rand_seed_det, $rand_seed_partition,
		   $nfold, $npartition,
		   $outdir_simobs, $outfile_head_simobs);
    printf("cmd = %s\n", $cmd);
    system($cmd);

    for(my $ipart = 1; $ipart <= $npartition; $ipart ++){
	my $logmu_lo = 0.0;
	my $logmu_up = 10.0;
	my $logmu_step = 1.0;
	my $beta_lo = 0.10;
	my $beta_up = 1.10;
	my $beta_step = 0.1;

	my ($mu_min, $beta_min);
	&RecEval($ipart, $respdir, $nfold,
		 $outdir_simobs, $outfile_head_simobs,
		 $logmu_lo, $logmu_up, $logmu_step,
		 $beta_lo, $beta_up, $beta_step,
		 \$mu_min, \$beta_min);
	printf("(mu_min, beta_min) = (%e, %e)\n", $mu_min, $beta_min);

    }
    printf("#=== %s: END ===#\n", $mipllib::sevar{'progname'});
} # main

sub Usage
{
    printf("usage: %s " .
	   " datalist\n",
	   $mipllib::sevar{'progname'});
    exit;
}

sub MkDataList
{
    if(4 != $#_){
        die;
    }
    my($datalist, $nfold, $outdir_simobs, $outfile_head_simobs, $ipart) = @_;
    open(DATALIST, ">$datalist");
    printf(DATALIST "# datafile  outdir\n");
    printf(DATALIST "%s/%s_sky_0001_obs_0001.fits  %s_sky_0001_obs_0001\n",
	   $outdir_simobs, $outfile_head_simobs, $outfile_head_simobs);
    printf(DATALIST "#\n");
    for(my $ifold = 0; $ifold < $nfold; $ifold ++){
	printf(DATALIST "%s/%s_sky_0001_obs_0001_part%2.2d_%2.2dfold%2.2d_tr.fits  " .
	       "%s_sky_0001_obs_0001_part%2.2d_%2.2dfold%2.2d_tr\n",
	       $outdir_simobs, $outfile_head_simobs, $ipart, $nfold, $ifold,
	       $outfile_head_simobs, $ipart, $nfold, $ifold);
    }
    close(DATALIST);
}

sub MkDataCvList
{
    if(5 != $#_){
        die;
    }
    my($data_cv_list, $nfold, $outdir_simobs, $outfile_head_simobs, $ipart, $outdir_emgist) = @_;
    open(DATALIST, ">$data_cv_list");
    printf(DATALIST "# fold_dir   recfile   valfile\n");
    for(my $ifold = 0; $ifold < $nfold; $ifold ++){
	printf(DATALIST "%s/%s_sky_0001_obs_0001_part%2.2d_%2.2dfold%2.2d_tr  " .
	       "rec_rec.fits  " .
	       "%s/%s_sky_0001_obs_0001_part%2.2d_%2.2dfold%2.2d_vl.fits\n",
	       $outdir_emgist, $outfile_head_simobs, $ipart, $nfold, $ifold,
	       $outdir_simobs, $outfile_head_simobs, $ipart, $nfold, $ifold);
    }
    close(DATALIST);
}

sub MkMulist
{
    if(3 != $#_){
        die;
    }
    my($mulist, $logmu_lo, $logmu_up, $logmu_step) = @_;
    open(MULIST, ">$mulist");
    my $logmu = $logmu_lo;
    while($logmu <= $logmu_up){
	printf(MULIST "%2.2e\n", pow(10.0, $logmu));
	$logmu += $logmu_step;
    }
    close(MULIST);
}

sub MkBetalist
{
    if(3 != $#_){
        die;
    }
    my($betalist, $beta_lo, $beta_up, $beta_step) = @_;
    open(BETALIST, ">$betalist");
    my $beta = $beta_lo;
    while($beta <= $beta_up){
	printf(BETALIST "%2.2e\n", $beta);
	$beta += $beta_step;
    }
    close(BETALIST);
}

sub RecEval
{
    if(12 != $#_){
        die;
    }
    my($ipart, $respdir, $nfold,
       $outdir_simobs, $outfile_head_simobs,
       $logmu_lo, $logmu_up, $logmu_step,
       $beta_lo, $beta_up, $beta_step,
       $mu_min_ref, $beta_min_ref) = @_;

    my $outdir_part = sprintf("part%2.2d", $ipart);
    my $cmd = sprintf("mkdir %s", $outdir_part);
    printf("cmd = %s\n", $cmd);
    system($cmd);
    $cmd = sprintf("mkdir %s/setup", $outdir_part);
    printf("cmd = %s\n", $cmd);
    system($cmd);

    my $datalist = sprintf("%s/setup/data.list", $outdir_part);
    &MkDataList($datalist, $nfold, $outdir_simobs, $outfile_head_simobs, $ipart);
    my $mulist = sprintf("%s/setup/mu.list", $outdir_part);
    my $betalist = sprintf("%s/setup/beta.list", $outdir_part);
    my $data_cv_list = sprintf("%s/setup/data_cv.list", $outdir_part);
    my $outdir_emgist = sprintf("%s/emgist_case2_lsig", $outdir_part);
    my $outdir_cv     = sprintf("%s/cv", $outdir_part);

    &MkMulist($mulist, $logmu_lo, $logmu_up, $logmu_step);
    &MkBetalist($betalist, $beta_lo, $beta_up, $beta_step);

    $cmd = sprintf("/home/morii/work/github/moriiism/srt/script/emgist_case2_lsig.pl  " .
		   "%s  %s  %s  %s  %s",
		   $datalist,
		   $mulist,
		   $betalist,
		   $respdir,
		   $outdir_emgist);
    printf("cmd = %s\n", $cmd);
    system($cmd);

    &MkDataCvList($data_cv_list, $nfold, $outdir_simobs, $outfile_head_simobs, $ipart, $outdir_emgist);

    #
    # evaluate
    #
    $cmd = sprintf("/home/morii/work/github/moriiism/srt/script/eval_cv.pl " .
		   "%s  %s  %s  %d  %s  %s" ,
		   $data_cv_list, $mulist, $betalist, $nfold, $outdir_cv, $respdir);
    printf("cmd = %s\n", $cmd);
    system($cmd);

    my $nbin_hist_logmu = int( ($logmu_up - $logmu_lo + 1.0e-10) / $logmu_step ) + 1;
    my $logmu_hist_lo = $logmu_lo - $logmu_step / 2.0;
    my $logmu_hist_up = $logmu_up + $logmu_step / 2.0;

    my $nbin_hist_beta = int( ($beta_up - $beta_lo + 1.0e-10) / $beta_step ) + 1;
    my $beta_hist_lo = $beta_lo - $beta_step / 2.0;
    my $beta_hist_up = $beta_up + $beta_step / 2.0;

    my $hist_info_file = sprintf("%s/setup/hist_info.dat", $outdir_part);
    open(HIST_INFO, ">$hist_info_file");
    printf(HIST_INFO "# nbin  lo  up  delta-bin  mode\n");
    printf(HIST_INFO "%d  %e     %e    none  none\n",
	   $nbin_hist_logmu, $logmu_hist_lo, $logmu_hist_up);
    printf(HIST_INFO "%d  %e     %e    none  none\n",
	   $nbin_hist_beta, $beta_hist_lo, $beta_hist_up);
    close(HIST_INFO);

    my $logmu_hist_index_lo = -0.5;
    my $logmu_hist_index_up = $nbin_hist_logmu - 0.5;
    my $beta_hist_index_lo = -0.5;
    my $beta_hist_index_up = $nbin_hist_beta - 0.5;

    my $hist_info_index_file = sprintf("%s/setup/hist_info_index.dat", $outdir_part);
    open(HIST_INFO, ">$hist_info_index_file");
    printf(HIST_INFO "# nbin  lo  up  delta-bin  mode\n");
    printf(HIST_INFO "%d  %e     %e    none  none\n",
	   $nbin_hist_logmu, $logmu_hist_index_lo, $logmu_hist_index_up);
    printf(HIST_INFO "%d  %e     %e    none  none\n",
	   $nbin_hist_beta, $beta_hist_index_lo, $beta_hist_index_up);
    close(HIST_INFO);

    ## find  zrange
    my $val_min = 1.0e10;
    my $val_max = -1.0e10;
    my $infile_smr = sprintf("%s/smr.txt", $outdir_cv);
    open(SMR, "<$infile_smr");
    while(my $line = <SMR>){
	chomp($line);
	if($line =~ /^\s*$/ or $line =~ /^\s*\#/ or $line =~ /^\s*\!/) {
            next;
        }
	my($imu, $ibeta, $mu, $beta, $val, $val_err) = split(' ', $line);
	if($val < $val_min){
	    $val_min = $val;
	}
	if($val > $val_max){
	    $val_max = $val;
	}
    }
    close(SMR);

    my $zrange_lo = $val_min - ($val_max - $val_min) * 0.1;
    my $zrange_up = $val_max + ($val_max - $val_min) * 0.1;
    my $outfile_head = "cvmap";

    $cmd = sprintf("/home/morii/work/github/moriiism/srt/cvmap/cvmap -- " .
		   "%s  %s  %s  %e  %e  %s  %s",
		   $infile_smr,
		   $hist_info_file,
		   $hist_info_index_file,
		   $zrange_lo,
		   $zrange_up,
		   $outdir_cv,
		   $outfile_head);
    printf("cmd = %s\n", $cmd);
    system($cmd);


    my $min_par_file = sprintf("%s/%s_out.dat", $outdir_cv, $outfile_head);
    open(MIN, "<$min_par_file");
    <MIN>;
    my($mu_min, $beta_min) = split(' ', <MIN>);
    close(MIN);

    $$mu_min_ref = $mu_min;
    $$beta_min_ref = $beta_min;
}
