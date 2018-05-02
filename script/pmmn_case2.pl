#!/usr/bin/perl

#
# pmmn_case2.pl
#
# 2018.05.02  M.Morii
#

use strict;
use warnings;
use POSIX;
no warnings qw(once);

$ENV{'LANG'} = "C";

my $mipllib = sprintf("%s/script/mipllib/mipllib.pl",
		      $ENV{'MITOOL'});
require($mipllib);
$mipllib::sevar{'progname'} = "pmmn_case2.pl";

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
		my $outfile_head = "rec";
		my $tol = 1.0e-10;
		my $nstep = 1000;
		my $lconst = 1.0;
		my $epsilon = 1.0e-15;
		$cmd = sprintf("/home/morii/work/github/moriiism/srt/pmmn/pmmn_case2  " .
			       "%s  %s  " . 
			       "%e  %e  " .
			       "%s  %s  " .
			       "%e  %d  %e  %e" ,
			       $respdir, $datafile, 
			       $mu_arr[$imu], $beta_arr[$ibeta], 
			       $outdir_this2, $outfile_head,
			       $tol, $nstep, $lconst, $epsilon);
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


