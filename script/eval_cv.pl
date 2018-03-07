#!/usr/bin/perl

#
# eval_cv.pl
#
# 2017.12.01  M.Morii
#

use strict;
use warnings;
use POSIX;
no warnings qw(once);

$ENV{'LANG'} = "C";

my $mipllib = sprintf("%s/script/mipllib/mipllib.pl",
		      $ENV{'MITOOL'});
require($mipllib);
$mipllib::sevar{'progname'} = "run_deconv.pl";

# main
{
    printf("#--- %s: START ---#\n", $mipllib::sevar{'progname'});

    my $cmd;
    # -- argument
    my $NARG = 3;
    my($mulist, $betalist, $nfold);
    printf("# of arg = %d\n", $#ARGV + 1);
    if( $#ARGV == $NARG - 1){
	my $iarg = 0;
	$mulist   = $ARGV[$iarg]; $iarg ++;
	$betalist = $ARGV[$iarg]; $iarg ++;
	$nfold     = $ARGV[$iarg]; $iarg ++;
    }else{
        printf("# of arg must be %d.\n", $NARG);
	&Usage();
    }
    my $prompt_arg = "arg";
    printf("%s: %s: mulist    = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $mulist);
    printf("%s: %s: betalist  = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $betalist);
    printf("%s: %s: nfold     = %d\n", $mipllib::sevar{'progname'}, $prompt_arg, $nfold);
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

    for(my $imu = 0; $imu <= $#mu_arr; $imu ++){
	for(my $ibeta = 0; $ibeta <= $#beta_arr; $ibeta ++){
	    my $outfile = sprintf("setup/cvfile_mu%.1e_beta%.1e.list",
				  $mu_arr[$imu], $beta_arr[$ibeta]);
	    open(OUT, ">$outfile");
	    for(my $ifold = 1; $ifold <= $nfold; $ifold ++){
		printf(OUT "out/double_star_off_029_031_cv_1e+03_tr_%dfold%d/mu%.1e/beta%.1e/deconv.fits  " .
		       "../data/20180228/doublestar_off_axis/doublestar_gimage_029_031_plus_gimage_029_029/cv_1e+03_vl_%dfold%d.img\n",
		       $nfold, $ifold, $mu_arr[$imu], $beta_arr[$ibeta], $nfold, $ifold);
	    }
	    close(OUT);
	}
    }

    printf("#=== %s: END ===#\n", $mipllib::sevar{'progname'});
} # main


sub Usage
{
    printf("usage: %s " .
	   " datalist  mulist  betalist  outdir\n",
	   $mipllib::sevar{'progname'});
    exit;
}


