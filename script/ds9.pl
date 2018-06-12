#!/usr/bin/perl

#
# ds9.pl
#
# 2018.05.14  M.Morii
#

use strict;
use warnings;
use POSIX;
no warnings qw(once);

$ENV{'LANG'} = "C";

my $mipllib = sprintf("%s/script/mipllib/mipllib.pl",
			$ENV{'MITOOL'});
require($mipllib);
$mipllib::sevar{'progname'} = "ds9.pl";

# main
{
    printf("#--- %s: START ---#\n", $mipllib::sevar{'progname'});

    my $cmd;
    # -- argument
    my $NARG = 6;
    my($datadir, $mulist, $betalist, $outdir, $scale_lo, $scale_up);
    printf("# of arg = %d\n", $#ARGV + 1);
    if( $#ARGV == $NARG - 1){
        my $iarg = 0;
        $datadir  = $ARGV[$iarg]; $iarg ++;
        $mulist   = $ARGV[$iarg]; $iarg ++;
        $betalist = $ARGV[$iarg]; $iarg ++;
        $outdir   = $ARGV[$iarg]; $iarg ++;
	$scale_lo = $ARGV[$iarg]; $iarg ++;
	$scale_up = $ARGV[$iarg]; $iarg ++;
    }else{
        printf("# of arg must be %d.\n", $NARG);
        &Usage();
    }
    my $prompt_arg = "arg";
    printf("%s: %s: datadir   = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $datadir);
    printf("%s: %s: mulist    = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $mulist);
    printf("%s: %s: betalist  = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $betalist);
    printf("%s: %s: outdir    = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $outdir);
    printf("%s: %s: scale_lo  = %e\n", $mipllib::sevar{'progname'}, $prompt_arg, $scale_lo);
    printf("%s: %s: scale_up  = %e\n", $mipllib::sevar{'progname'}, $prompt_arg, $scale_up);
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
    
    my $filelist_str = "";

    for(my $imu = 0; $imu <= $#mu_arr; $imu ++){
	## for(my $ibeta = 0; $ibeta <= $#beta_arr; $ibeta ++){
	for(my $ibeta = $#beta_arr; $ibeta >= 0; $ibeta --){
	    
	    my $file_this = sprintf("%s/mu%1.1e/beta%1.2e/rec_rec.fits",
				    $datadir, 
				    $mu_arr[$imu], $beta_arr[$ibeta]);
	    $filelist_str .= " " . $file_this;
	}
	$filelist_str .= " ";
    }
    printf("%s\n", $filelist_str);

    # -cmap value <contrast> <bias> 

    my $outimg = sprintf("%s/rec_all.png", $outdir);

    $cmd = sprintf("ds9 -xpa no " .
		   "-geometry 1400x950 " .
		   "-tile grid layout %d %d " .
		   "-zoom to fit " .
		   "-scale mode minmax " .
		   "-scale scope global " .
		   "-scale limits %e %e " .
		   "-cmap value 5.0 0.1 " .
		   "%s " .
		   "-saveimage png %s " ,
		   $#beta_arr + 1, $#mu_arr + 1,
		   $scale_lo, $scale_up,
		   $filelist_str,
		   $outimg);
    printf("cmd = %s\n", $cmd);
    system($cmd);

    printf("#=== %s: END ===#\n", $mipllib::sevar{'progname'});
} # main


sub Usage
{
    printf("usage: %s " .
           " datalist  mulist  betalist  outdir\n",
           $mipllib::sevar{'progname'});
    exit;
}
