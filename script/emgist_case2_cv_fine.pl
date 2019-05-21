#!/usr/bin/perl

#
# emgist_case2_cv_fine.pl
#
# 2018.08.20  M.Morii
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
$mipllib::sevar{'progname'} = "emgist_case2_cv_fine.pl";

# main
{
    printf("#--- %s: START ---#\n", $mipllib::sevar{'progname'});

    my $cmd;
    # -- argument
    my $NARG = 5;
    my($datalist, $data_cv_fine_list, $respdir, $mu_init, $beta_init);
    printf("# of arg = %d\n", $#ARGV + 1);
    if( $#ARGV == $NARG - 1){
	my $iarg = 0;
	$datalist     = $ARGV[$iarg]; $iarg ++;
	$data_cv_fine_list = $ARGV[$iarg]; $iarg ++;
	$respdir      = $ARGV[$iarg]; $iarg ++;
	$mu_init      = $ARGV[$iarg]; $iarg ++;
	$beta_init    = $ARGV[$iarg]; $iarg ++;
    }else{
        printf("# of arg must be %d.\n", $NARG);
	&Usage();
    }
    my $prompt_arg = "arg";
    printf("%s: %s: datalist      = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $datalist);
    printf("%s: %s: data_cv_fine_list  = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $data_cv_fine_list);
    printf("%s: %s: respdir       = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $respdir);
    printf("%s: %s: mu_init       = %e\n", $mipllib::sevar{'progname'}, $prompt_arg, $mu_init);
    printf("%s: %s: beta_init     = %e\n", $mipllib::sevar{'progname'}, $prompt_arg, $beta_init);
    # == argument

    my $logmu_init = log10($mu_init);
    my $logmu_step = 1.0;
    my $beta_step = 0.1;
    my $mulist_fine = "setup/mu_fine2.list";
    my $logmu_lo_fine = $logmu_init - $logmu_step;
    my $logmu_up_fine = $logmu_init + $logmu_step;
    my $logmu_step_fine = $logmu_step / 5.0;
    my $betalist_fine = "setup/beta_fine2.list";
    my $beta_lo_fine = $beta_init - $beta_step;
    my $beta_up_fine = $beta_init + $beta_step;
    my $beta_step_fine = $beta_step / 5.0;

    my $mu_min_fine = 0.0;
    my $beta_min_fine = 0.0;
    my $outdir_emgist_fine = "emgist_case2_lsig_fine2";
    my $outdir_cv_fine     = "cv_fine2";
    &RecEval($datalist, $data_cv_fine_list, $respdir,
	     $mulist_fine, $logmu_lo_fine, $logmu_up_fine, $logmu_step_fine,
	     $betalist_fine, $beta_lo_fine, $beta_up_fine, $beta_step_fine,
	     $outdir_emgist_fine, $outdir_cv_fine,
	     \$mu_min_fine, \$beta_min_fine);

    printf("(mu_min_fine, beta_min_fine) = (%e, %e)\n", $mu_min_fine, $beta_min_fine);

    printf("#=== %s: END ===#\n", $mipllib::sevar{'progname'});
} # main

sub Usage
{
    printf("usage: %s " .
	   " datalist\n",
	   $mipllib::sevar{'progname'});
    exit;
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
    if(14 != $#_){
        die;
    }
    my($datalist, $data_cv_list, $respdir,
       $mulist, $logmu_lo, $logmu_up, $logmu_step,
       $betalist, $beta_lo, $beta_up, $beta_step,
       $outdir_emgist, $outdir_cv,
       $mu_min_ref, $beta_min_ref) = @_;

    &MkMulist($mulist, $logmu_lo, $logmu_up, $logmu_step);
    &MkBetalist($betalist, $beta_lo, $beta_up, $beta_step);

    my $cmd = sprintf("/home/morii/work/github/moriiism/srt/script/emgist_case2_lsig.pl  " .
		      "%s  %s  %s  %s  %s",
		      $datalist,
		      $mulist,
		      $betalist,
		      $respdir,
		      $outdir_emgist);
    printf("cmd = %s\n", $cmd);
    system($cmd);

    #
    # evaluate
    #
    my $nfold = 5;
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

    my $hist_info_file = "setup/hist_info.dat";
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

    my $hist_info_index_file = "setup/hist_info_index.dat";
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
