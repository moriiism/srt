#!/usr/bin/perl

#
# smrmap_pre.pl
#
# 2018.06.08  M.Morii
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
$mipllib::sevar{'progname'} = "smrmap_pre.pl";

# main
{
    printf("#--- %s: START ---#\n", $mipllib::sevar{'progname'});

    my $cmd;
    # -- argument
    my $NARG = 4;
    my($datadir, $mulist, $betalist, $outdir);
    printf("# of arg = %d\n", $#ARGV + 1);
    if( $#ARGV == $NARG - 1){
	my $iarg = 0;
	$datadir  = $ARGV[$iarg]; $iarg ++;
	$mulist   = $ARGV[$iarg]; $iarg ++;
	$betalist = $ARGV[$iarg]; $iarg ++;
	$outdir   = $ARGV[$iarg]; $iarg ++;
    }else{
        printf("# of arg must be %d.\n", $NARG);
	&Usage();
    }
    my $prompt_arg = "arg";
    printf("%s: %s: datadir   = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $datadir);
    printf("%s: %s: mulist    = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $mulist);
    printf("%s: %s: betalist  = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $betalist);
    printf("%s: %s: outdir    = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $outdir);
    # == argument

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

    my $outfile = sprintf("%s/rec_last_all.dat", $outdir);
    open(OUT, ">$outfile");
    for(my $imu = 0; $imu <= $#mu_arr; $imu ++){
	for(my $ibeta = 0; $ibeta <= $#beta_arr; $ibeta ++){
	    my $infile = sprintf("%s/mu%1.1e/beta%1.2e/rec_last.dat",
				 $datadir,
				 $mu_arr[$imu], $beta_arr[$ibeta]);
	    $cmd = sprintf("cat %s", $infile);
	    my $cmdout = `$cmd`;
	    chomp($cmdout);
	    printf(OUT "%d  %d  %s\n", $imu, $ibeta, $cmdout);
	}
    }
    close(OUT);

    printf("#=== %s: END ===#\n", $mipllib::sevar{'progname'});
} # main


sub Usage
{
    printf("usage: %s " .
	   " datadir  mulist  betalist  outdir\n",
	   $mipllib::sevar{'progname'});
    exit;
}


