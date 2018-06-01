#!/usr/bin/perl

#
# emgist_case2.pl
#
# 2018.05.31  M.Morii
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
$mipllib::sevar{'progname'} = "emgist_case2.pl";

# main
{
    printf("#--- %s: START ---#\n", $mipllib::sevar{'progname'});

    my $cmd;
    # -- argument
    my $NARG = 4;
    my($datalist, $mulist, $betalist, $outdir);
    printf("# of arg = %d\n", $#ARGV + 1);
    if( $#ARGV == $NARG - 1){
	my $iarg = 0;
	$datalist = $ARGV[$iarg]; $iarg ++;
	$mulist   = $ARGV[$iarg]; $iarg ++;
	$betalist = $ARGV[$iarg]; $iarg ++;
	$outdir   = $ARGV[$iarg]; $iarg ++;
    }else{
        printf("# of arg must be %d.\n", $NARG);
	&Usage();
    }
    my $prompt_arg = "arg";
    printf("%s: %s: datalist  = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $datalist);
    printf("%s: %s: mulist    = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $mulist);
    printf("%s: %s: betalist  = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $betalist);
    printf("%s: %s: outdir    = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $outdir);
    # == argument

    my $respdir = "/home/morii/work/maeda/data/20170428b/model";

    my @mu_arr = ();
    my $iline = 0;
    open(LIST, "<$mulist");
    while(my $line = <LIST>){
	chomp($line);
        if($line =~ /^\s*$/ or $line =~ /^\s*\#/ or $line =~ /^\s*\!/) {
            next;
        }
	$mu_arr[$iline] = $line; $iline ++;
    }
    close(LIST);

    my @beta_arr = ();
    $iline = 0;
    open(LIST, "<$betalist");
    while(my $line = <LIST>){
	chomp($line);
        if($line =~ /^\s*$/ or $line =~ /^\s*\#/ or $line =~ /^\s*\!/) {
            next;
        }
	$beta_arr[$iline] = $line; $iline ++;
    }
    close(LIST);

    open(LIST, "<$datalist");
    while(my $line = <LIST>){
	chomp($line);
        if($line =~ /^\s*$/ or $line =~ /^\s*\#/ or $line =~ /^\s*\!/) {
            next;
        }
	my ($datafile, $outdir_this) = split(' ', $line);

	for(my $imu = 0; $imu <= $#mu_arr; $imu ++){
	    for(my $ibeta = 0; $ibeta <= $#beta_arr; $ibeta ++){

		my $outdir_this2 = sprintf("%s/%s/mu%1.1e/beta%1.2e",
					   $outdir, $outdir_this,
					   $mu_arr[$imu], $beta_arr[$ibeta]);
		if(! -e $outdir_this2){
		    $cmd = sprintf("mkdir -p %s", $outdir_this2);
		    printf("cmd = %s\n", $cmd);
		    system($cmd);
		}
		my $logfile = sprintf("%s/emgist_case2.log", $outdir_this2);

		my $skyfile = "";
		if(0 == $ibeta){
		    $skyfile = "none";
		} else {
		    $skyfile = sprintf("%s/%s/mu%1.1e/beta%1.2e/rec_rec.fits",
				       $outdir, $outdir_this,
				       $mu_arr[$imu], $beta_arr[$ibeta - 1]);
		}

		### temp
		$skyfile = "none";

		## check flag file
		my $flag_file_openblas = "setup/openblas_on";
		my $exeprog = sprintf("emgist");
		if(-e $flag_file_openblas){
		    $exeprog = sprintf("emgist_openblas");
		}

		my $flag_line_search = 1;

		my $outfile_head = "rec";
		my $nem    = 3000;
		my $tol_em = 1.0e-10;
		my $npm    = 1000;
		my $tol_pm = 1.0e-10;
		my $lconst = 1.0e-5;
		my $epsilon = 1.0e-15;
		$cmd = sprintf("stdbuf -oL -eL /home/morii/work/github/moriiism/srt/emgist_case2/%s  " .
			       "%s  %s  %s  " . 
			       "%e  %e  " .
			       "%s  %s  " .
			       "%d  %e  %d  %e  %d  %e  %e " .
			       "> %s 2>&1" ,
			       $exeprog,
			       $respdir, $datafile, $skyfile,
			       $mu_arr[$imu], $beta_arr[$ibeta], 
			       $outdir_this2, $outfile_head,
			       $nem, $tol_em, $npm, $tol_pm, $flag_line_search, $lconst, $epsilon,
			       $logfile);
		printf("cmd = %s\n", $cmd);
		system($cmd);
	    }
	}
    }
    close(LIST);

    printf("#=== %s: END ===#\n", $mipllib::sevar{'progname'});
} # main


sub Usage
{
    printf("usage: %s " .
	   " datalist  mulist  betalist  outdir\n",
	   $mipllib::sevar{'progname'});
    exit;
}


