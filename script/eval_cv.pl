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
$mipllib::sevar{'progname'} = "eval_cv.pl";

# main
{
    printf("#--- %s: START ---#\n", $mipllib::sevar{'progname'});

    my $cmd;
    # -- argument
    my $NARG = 6;
    my($cvlist, $mulist, $betalist, $nfold, $outdir, $respdir);
    printf("# of arg = %d\n", $#ARGV + 1);
    if( $#ARGV == $NARG - 1){
	my $iarg = 0;
	$cvlist   = $ARGV[$iarg]; $iarg ++;
	$mulist   = $ARGV[$iarg]; $iarg ++;
	$betalist = $ARGV[$iarg]; $iarg ++;
	$nfold    = $ARGV[$iarg]; $iarg ++;
	$outdir   = $ARGV[$iarg]; $iarg ++;
	$respdir  = $ARGV[$iarg]; $iarg ++;
    }else{
        printf("# of arg must be %d.\n", $NARG);
	&Usage();
    }
    my $prompt_arg = "arg";
    printf("%s: %s: cvlist    = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $cvlist);
    printf("%s: %s: mulist    = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $mulist);
    printf("%s: %s: betalist  = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $betalist);
    printf("%s: %s: nfold     = %d\n", $mipllib::sevar{'progname'}, $prompt_arg, $nfold);
    printf("%s: %s: outdir    = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $outdir);
    printf("%s: %s: respdir   = %s\n", $mipllib::sevar{'progname'}, $prompt_arg, $respdir);
    # == argument

    my $outdir_setup = sprintf("%s/setup", $outdir);
    if(! -e $outdir_setup){
	$cmd = sprintf("mkdir -p %s", $outdir_setup);
	printf("cmd = %s\n", $cmd);
	system($cmd);
    }

    my @fold_dir_arr = ();
    my @recfile_arr  = ();
    my @valfile_arr  = ();
    my $iline = 0;
    open(LIST, "<$cvlist");
    while(my $line = <LIST>){
	chomp($line);
        if($line =~ /^\s*$/ or $line =~ /^\s*\#/ or $line =~ /^\s*\!/) {
            next;
        }
	($fold_dir_arr[$iline], $recfile_arr[$iline], $valfile_arr[$iline]) = split(' ', $line);
	$iline ++;
    }
    close(LIST);

    my @mu_arr = ();
    $iline = 0;
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

    my $smrfile = sprintf("%s/smr.txt", $outdir);
    open(SMR, ">$smrfile");

    for(my $imu = 0; $imu <= $#mu_arr; $imu ++){
	for(my $ibeta = 0; $ibeta <= $#beta_arr; $ibeta ++){
	    my $filelist = sprintf("%s/mu%.2e_beta%.2e.list",
				  $outdir_setup, $mu_arr[$imu], $beta_arr[$ibeta]);
	    open(OUT, ">$filelist");
	    for(my $ifold = 0; $ifold < $nfold; $ifold ++){
		printf(OUT "%s/mu%.2e/beta%.2e/%s  %s\n" ,
		       $fold_dir_arr[$ifold], $mu_arr[$imu], $beta_arr[$ibeta], $recfile_arr[$ifold],
		       $valfile_arr[$ifold]);
	    }
	    close(OUT);

	    my $outfile_head = sprintf("mu%.2e_beta%.2e",
				       $mu_arr[$imu], $beta_arr[$ibeta]);
	    $cmd = sprintf("/home/morii/work/github/moriiism/srt/eval/eval_cv  " .
			   "%s  %s  %s  %s",
			   $respdir, $filelist,
			   $outdir, $outfile_head);
	    printf("cmd = %s\n", $cmd);
	    system($cmd);
	    
	    my $heldistfile = sprintf("%s/%s_heldist.txt", $outdir, $outfile_head);
	    $cmd = sprintf("cat %s", $heldistfile);
	    my $heldist = `$cmd`;
	    chomp($heldist);
	    printf(SMR "%d  %d  %e  %e  %s\n", 
		   $imu, $ibeta,
		   $mu_arr[$imu], $beta_arr[$ibeta], $heldist);
	}
    }
    close(SMR);

    printf("#=== %s: END ===#\n", $mipllib::sevar{'progname'});
} # main


sub Usage
{
    printf("usage: %s " .
	   " cvlist  mulist  betalist  nfold\n",
	   $mipllib::sevar{'progname'});
    exit;
}

