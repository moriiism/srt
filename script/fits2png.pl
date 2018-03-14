#!/usr/bin/perl

#
# fits2png.pl
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
$mipllib::sevar{'progname'} = "fits2png.pl";

# main
{
    printf("#--- %s: START ---#\n", $mipllib::sevar{'progname'});

    my $cmd;
    # -- argument
    my $NARG = 1;
    my($datalist);
    printf("# of arg = %d\n", $#ARGV + 1);
    if( $#ARGV == $NARG - 1){
	my $iarg = 0;
	$datalist = $ARGV[$iarg]; $iarg ++;
    }else{
        printf("# of arg must be %d.\n", $NARG);
	&Usage();
    }
    my $prompt_arg = "arg";
    printf("%s: %s: datalist  = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $datalist);
    # == argument

    open(LIST, "<$datalist");
    while(my $line = <LIST>){
	chomp($line);
        if($line =~ /^\s*$/ or $line =~ /^\s*\#/ or $line =~ /^\s*\!/) {
            next;
        }
	my $infile = $line;
	$cmd = sprintf("dirname %s", $infile);
	my $outdir = `$cmd`;
	chomp($outdir);

	$cmd = sprintf("basename %s", $infile);
	my $outfile_head = `$cmd`;
	chomp($outfile_head);

	my $outpng = sprintf("%s/%s.png", $outdir, $outfile_head);
	$cmd = sprintf("ds9  %s  " .
		       "-geometry 640x900 " .
		       "-zoom to fit " .
		       "-scale mode minmax " .
		       "-colorbar no " .
		       "-saveimage png " .
		       "%s " . 
		       "-quit" ,
		       $infile, $outpng);
	printf("cmd = %s\n", $cmd);
	system($cmd);
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
