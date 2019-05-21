#!/usr/bin/perl

#
# emgist_case2_sim_npeak.pl
#
# 2018.08.22  M.Morii
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
$mipllib::sevar{'progname'} = "emgist_case2_sim_npeak.pl";

# main
{
    printf("#--- %s: START ---#\n", $mipllib::sevar{'progname'});

    my $cmd;
    # -- argument
    my $NARG = 6;
    my($infile, $nevt, $nobs, $respdir, $mu, $beta);

    printf("# of arg = %d\n", $#ARGV + 1);
    if( $#ARGV == $NARG - 1){
	my $iarg = 0;
	$infile       = $ARGV[$iarg]; $iarg ++;
	$nevt         = $ARGV[$iarg]; $iarg ++;
	$nobs         = $ARGV[$iarg]; $iarg ++;
	$respdir      = $ARGV[$iarg]; $iarg ++;
	$mu           = $ARGV[$iarg]; $iarg ++;
	$beta         = $ARGV[$iarg]; $iarg ++;
    }else{
        printf("# of arg must be %d.\n", $NARG);
	&Usage();
    }
    my $prompt_arg = "arg";
    printf("%s: %s: infile        = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $infile);
    printf("%s: %s: nevt          = %d\n", $mipllib::sevar{'progname'}, $prompt_arg, $nevt);
    printf("%s: %s: nobs          = %d\n", $mipllib::sevar{'progname'}, $prompt_arg, $nobs);
    printf("%s: %s: respdir       = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $respdir);
    printf("%s: %s: mu            = %e\n", $mipllib::sevar{'progname'}, $prompt_arg, $mu);
    printf("%s: %s: beta          = %e\n", $mipllib::sevar{'progname'}, $prompt_arg, $beta);
    # == argument

    my $outdir_simobs = "simobs";
    my $outfile_head_simobs = "simobs";

    for(my $iobs = 0; $iobs < $nobs; $iobs ++){
	# simulation
	my $rand_seed_sky = 1;
	my $rand_seed_det = 1 + $iobs;
	my $rand_seed_partition = 1;
	my $nfold = 1;
	my $npartition = 1;
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
    }

    my $npeak_min = 1000;
    my $npeak_max = -1000;
    my @num_npeak_arr = ();
    for(my $inum = 0; $inum < 100; $inum ++){
	$num_npeak_arr[$inum] = 0;
    }
    for(my $iobs = 1; $iobs <= $nobs; $iobs ++){
	my $npeak = 0;
	&RecNpeak($iobs, $mu, $beta, $respdir, $outdir_simobs, $outfile_head_simobs, \$npeak);
	printf("npeak = %d\n", $npeak);
	if($npeak < $npeak_min){
	    $npeak_min = $npeak;
	}
	if($npeak > $npeak_max){
	    $npeak_max = $npeak;
	}
	$num_npeak_arr[$npeak] ++;
    }
    printf("npeak_min, npeak_max = %d, %d\n", $npeak_min, $npeak_max);
    printf("num_npeak_arr[0] = %d\n", $num_npeak_arr[0]);
    printf("num_npeak_arr[1] = %d\n", $num_npeak_arr[1]);
    printf("num_npeak_arr[2] = %d\n", $num_npeak_arr[2]);
    printf("num_npeak_arr[3] = %d\n", $num_npeak_arr[3]);
    my $rate_success = $num_npeak_arr[2] / $nobs;
    printf("success rate = %e\n", $rate_success);

    $cmd = "mkdir sim_npeak";
    printf("cmd = %s\n", $cmd);
    system($cmd);
    my $outdat = sprintf("sim_npeak/npeak.dat");
    open(OUT, ">$outdat");
    printf(OUT "%e\n", $rate_success);
    close(OUT);

    printf("#=== %s: END ===#\n", $mipllib::sevar{'progname'});
} # main

sub Usage
{
    printf("usage: %s " .
	   " datalist\n",
	   $mipllib::sevar{'progname'});
    exit;
}

sub RecNpeak
{
    if(6 != $#_){
        die;
    }
    my $cmd;
    my($iobs, $mu, $beta, $respdir,
       $outdir_simobs, $outfile_head_simobs,
       $npeak_ref) = @_;
    
    my $outdir = sprintf("sim_npeak%3.3d", $iobs);
    if(! -e $outdir){
	$cmd = sprintf("mkdir %s", $outdir);
	printf("cmd = %s\n", $cmd);
	system($cmd);
    }

    my $datafile = sprintf("%s/%s_sky_0001_obs_%4.4d.fits",
			   $outdir_simobs, $outfile_head_simobs, $iobs);
    my $logfile = sprintf("%s/emgist_case2_lsig.log", $outdir);

    my $skyfile = "none";
    my $exeprog = sprintf("emgist");
    my $flag_line_search = 1;

    my $outfile_head = "rec";
    my $nem    = 3000;
    my $tol_em = 1.0e-10;
    my $npm    = 1000;
    ## my $tol_pm = 1.0e-10;
    my $tol_pm = 1.0e-5;
    my $tol_diff_l_var = 5.0e-2;
    my $lconst = 1.0e-5;
    my $epsilon = 1.0e-15;
    $cmd = sprintf("stdbuf -oL -eL /home/morii/work/github/moriiism/srt/emgist_case2_lsig/%s  " .
		   "%s  %s  %s  " . 
		   "%e  %e  " .
		   "%s  %s  " .
		   "%d  %e  %d  %e  %e  %d  %e  %e " .
		   "> %s 2>&1" ,
		   $exeprog,
		   $respdir, $datafile, $skyfile,
		   $mu, $beta, 
		   $outdir, $outfile_head,
		   $nem, $tol_em, $npm, $tol_pm, $tol_diff_l_var, $flag_line_search, $lconst, $epsilon,
		   $logfile);
    printf("cmd = %s\n", $cmd);
    system($cmd);

    my $recfile = sprintf("%s/rec_rec.fits", $outdir);
    my $xpos = 45;
    my $xpos = 30;
    my $ypos = 30;
    my $theta = 90.0;
    my $y_lo = -5;
    my $y_up = 5;
    my $nbinx_new = 70;
    my $xlo_new = -35.0;
    my $xup_new = 35.0;
    my $significance = 3.0;
    my $outdir_npeak = sprintf("%s", $outdir);
    $outfile_head = "npeak";

    $cmd = sprintf("/home/morii/work/github/moriiism/srt/npeak/npeak " .
		   "-- " .
		   "%s  %e  %e  %e " .
		   "%e  %e  %d  %e  %e " .
		   "%e  %s  %s  ",
		   $recfile, $xpos, $ypos, $theta,
		   $y_lo, $y_up, $nbinx_new, $xlo_new, $xup_new,
		   $significance, $outdir, $outfile_head);
    printf("cmd = %s\n", $cmd);
    system($cmd);

    my $npeak_par_file = sprintf("%s/%s_out.dat", $outdir, $outfile_head);
    open(NPEAK, "<$npeak_par_file");
    my $npeak = <NPEAK>;
    chomp($npeak);
    close(NPEAK);
    $$npeak_ref = $npeak;
}
