#!/usr/bin/perl

#
# mkallfits.pl
#
# 2018.06.18  M.Morii
#

use strict;
use warnings;
use POSIX;
no warnings qw(once);

$ENV{'LANG'} = "C";

my $mipllib = sprintf("%s/script/mipllib/mipllib.pl",
		      $ENV{'MITOOL'});
require($mipllib);
$mipllib::sevar{'progname'} = "mkallfits.pl";

# main
{
    printf("#--- %s: START ---#\n", $mipllib::sevar{'progname'});

    my $cmd;
    # -- argument
    my $NARG = 4;
    my($indata, $nevt, $nout, $outdir);
    printf("# of arg = %d\n", $#ARGV + 1);
    if( $#ARGV == $NARG - 1){
	my $iarg = 0;
	$indata = $ARGV[$iarg]; $iarg ++;
	$nevt   = $ARGV[$iarg]; $iarg ++;
	$nout   = $ARGV[$iarg]; $iarg ++;
	$outdir = $ARGV[$iarg]; $iarg ++;
    }else{
        printf("# of arg must be %d.\n", $NARG);
	&Usage();
    }
    my $prompt_arg = "arg";
    printf("%s: %s: indata  = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $indata);
    printf("%s: %s: nevt    = %d\n", $mipllib::sevar{'progname'}, $prompt_arg, $nevt);
    printf("%s: %s: nout    = %d\n", $mipllib::sevar{'progname'}, $prompt_arg, $nout);
    printf("%s: %s: outdir  = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $outdir);
    # == argument

    if(! -e $outdir){
	$cmd = sprintf("mkdir -p %s", $outdir);
	printf("cmd = %s\n", $cmd);
	system($cmd);
    }

    for(my $iout = 0; $iout < $nout; $iout ++){

	my $outfile_tmp = sprintf("%s/all_%.0e_%2.2d.fits", 
				  $outdir, $nevt, $iout);
	my $outfile     = sprintf("%s/all_%.0e_%2.2d.img", 
				  $outdir, $nevt, $iout);

	my $irow_st = $iout * $nevt + 1;
	my $irow_ed = $irow_st + $nevt;
	$cmd = sprintf("ftselect infile=%s outfile=%s " .
		       "expression=\"%d <= #row && #row<= %d \"",
		       $indata, $outfile_tmp,
		       $irow_st, $irow_ed);
	printf("cmd = %s\n", $cmd);
	system($cmd);

	$cmd = sprintf("f2dhisto clobber=yes  %s %s " .
		       "1 1 sudarexpos sudareypos \"-30,30\" \"-30,30\"",
		       $outfile_tmp, $outfile);
	printf("cmd = %s\n", $cmd);
	system($cmd);
    }

    close(LIST);

    printf("#=== %s: END ===#\n", $mipllib::sevar{'progname'});
} # main

sub Usage
{
    printf("usage: %s " .
	   "indata  nevt  nout  outdir\n",
	   $mipllib::sevar{'progname'});
    exit;
}
