#!/usr/bin/perl

#
# mkcvfits.pl
#
# 2018.03.02  M.Morii
#

use strict;
use warnings;
use POSIX;
no warnings qw(once);

$ENV{'LANG'} = "C";

my $mipllib = sprintf("%s/script/mipllib/mipllib.pl",
		      $ENV{'MITOOL'});
require($mipllib);
$mipllib::sevar{'progname'} = "mkcvfits.pl";

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

    my @nevt_arr = (1.0e3, 1.0e4);
    my $nfold = 5;

    my $iline = 0;
    open(LIST, "<$datalist");
    while(my $line = <LIST>){
	chomp($line);
        if($line =~ /^\s*$/ or $line =~ /^\s*\#/ or $line =~ /^\s*\!/) {
            next;
        }
	my $infile = $line;
	my $dir = dirname($infile);
	printf("dir = %s\n", $dir);

	# all data
	for(my $ievt_index = 0; $ievt_index <= $#nevt_arr; $ievt_index ++){
	    my $outfile_tmp = sprintf("%s/all_%.0e.fits", 
				      $dir, $nevt_arr[$ievt_index]);
	    my $outfile     = sprintf("%s/all_%.0e.img", 
				      $dir, $nevt_arr[$ievt_index]);

	    my $irow_st = 1;
	    my $irow_ed = $nevt_arr[$ievt_index];
	    $cmd = sprintf("ftselect infile=%s outfile=%s " .
			   "expression=\"%d <= #row && #row<= %d \"",
			   $infile, $outfile_tmp,
			   $irow_st, $irow_ed);
	    printf("cmd = %s\n", $cmd);
	    system($cmd);

	    $cmd = sprintf("f2dhisto clobber=yes  %s %s " .
			   "1 1 sudarexpos sudareypos \"-30,30\" \"-30,30\"",
			   $outfile_tmp, $outfile);
	    printf("cmd = %s\n", $cmd);
	    system($cmd);
	} 

	# cv data
	for(my $ievt_index = 0; $ievt_index <= $#nevt_arr; $ievt_index ++){
	    for(my $ifold = 1; $ifold <= $nfold; $ifold ++){
		my $infile_cv   = sprintf("%s/all_%.0e.fits",
					  $dir, $nevt_arr[$ievt_index]);
		my $outfile_tr_tmp = sprintf("%s/cv_%.0e_tr_%dfold%d.fits",
					     $dir, $nevt_arr[$ievt_index], $nfold, $ifold);
		my $outfile_tr     = sprintf("%s/cv_%.0e_tr_%dfold%d.img", 
					     $dir, $nevt_arr[$ievt_index], $nfold, $ifold);

		my $outfile_vl_tmp = sprintf("%s/cv_%.0e_vl_%dfold%d.fits",
					     $dir, $nevt_arr[$ievt_index], $nfold, $ifold);
		my $outfile_vl     = sprintf("%s/cv_%.0e_vl_%dfold%d.img", 
					     $dir, $nevt_arr[$ievt_index], $nfold, $ifold);

		my $irow_vl_st = $nevt_arr[$ievt_index] / $nfold * ($ifold - 1) + 1;
		my $irow_vl_ed = $nevt_arr[$ievt_index] / $nfold * $ifold;
		my $irow_tr_rm_st = $irow_vl_st;
		my $irow_tr_rm_ed = $irow_vl_ed;

		$cmd = sprintf("ftselect infile=%s outfile=%s " .
			       "expression=\"!(%d <= #row && #row<= %d)\"",
			       $infile_cv, $outfile_tr_tmp,
			       $irow_tr_rm_st, $irow_tr_rm_ed);
		printf("cmd = %s\n", $cmd);
		system($cmd);
		$cmd = sprintf("f2dhisto clobber=yes  %s %s " .
			       "1 1 sudarexpos sudareypos \"-30,30\" \"-30,30\"",
			       $outfile_tr_tmp, $outfile_tr);
		printf("cmd = %s\n", $cmd);
		system($cmd);

		$cmd = sprintf("ftselect infile=%s outfile=%s " .
			       "expression=\"%d <= #row && #row<= %d \"",
			       $infile_cv, $outfile_vl_tmp,
			       $irow_vl_st, $irow_vl_ed);
		printf("cmd = %s\n", $cmd);
		system($cmd);
		$cmd = sprintf("f2dhisto clobber=yes  %s %s " .
			       "1 1 sudarexpos sudareypos \"-30,30\" \"-30,30\"",
			       $outfile_vl_tmp, $outfile_vl);
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
	   " datalist\n",
	   $mipllib::sevar{'progname'});
    exit;
}
