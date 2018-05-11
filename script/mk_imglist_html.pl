#!/usr/bin/perl

#
# mk_imglist_html.pl
# 
# make image list by html
#

use strict;
use warnings;
use Cwd;
use File::Basename;
use File::Temp;
use File::Spec;
use POSIX;
use Getopt::Long;
no warnings qw(once);

$ENV{'LANG'} = "C";

my $mipllib = sprintf("%s/script/mipllib/mipllib.pl",
		      $ENV{'MITOOL'});
require($mipllib);
$mipllib::sevar{'progname'} = "mk_imglist_html.pl";

# main
{
    my $cmd;

    # -- argument
    my $NARG = 4;
    my($datalist, $mulist, $betalist, $outdir);  # relpath from $procdir or fullpath
    printf("# of arg = %d\n", $#ARGV + 1);
    if( $#ARGV == $NARG - 1){
        my $iarg = 0;
        $datalist   = $ARGV[$iarg]; $iarg ++;
        $mulist     = $ARGV[$iarg]; $iarg ++;
	$betalist   = $ARGV[$iarg]; $iarg ++;
	$outdir     = $ARGV[$iarg]; $iarg ++;
    }else{
        printf("# of arg must be %d.\n", $NARG);
        &Usage();
    }
    my $prompt_arg = "arg";
    printf("%s: %s: datalist  = %s\n", 
	   $mipllib::sevar{'progname'}, $prompt_arg, $datalist);
    printf("%s: %s: mulist    = %s\n", 
	   $mipllib::sevar{'progname'}, $prompt_arg, $mulist);
    printf("%s: %s: betalist  = %s\n", 
	   $mipllib::sevar{'progname'}, $prompt_arg, $betalist);
    printf("%s: %s: outdir    = %s\n", 
	   $mipllib::sevar{'progname'}, $prompt_arg, $outdir);
    # == argument

    chdir($mipllib::sevar{'procdir'});
    my $procdir_full = cwd;
    my $imglist_full = File::Spec->rel2abs($imglist, $procdir_full);
    my $outhtml_full = File::Spec->rel2abs($outhtml, $procdir_full);

    printf("imglist_full = %s\n", $imglist_full);
    printf("outhtml_full = %s\n", $outhtml_full);

    &MkImageListHtml($imglist_full, $outhtml_full, "temp");

} # main


sub Usage
{
    my $func_name = "Usage";
    if(-1 != $#_){
        printf("%s: %s: error.\n", $mipllib::sevar{'progname'}, $func_name);
        die;
    }
    printf("usage: %s [--help=(0)] [--verbose=(0)] [--debug=(0)] " .
           "[--outflag=(overwrite)] [--procdir=(.)] " .
           "imglist  outhtml\n", $mipllib::sevar{'progname'});
    exit;
}



sub MkImageListHtml
{
    if(2 != $#_){
        die;
    }
    my($imglist_full, $outhtml_full, $caption) = @_;
    my $cmd;

    # load imagelist
    my @list = ();
    open(IN, "<$imglist_full") or die;
    while(my $line = <IN>){
        chomp($line);
        if($line =~ /^\s*$/ or $line =~ /^\s*\#/ or $line =~ /^\s*\!/) {
            next;
        }
	push(@list, $line);
    }
    close(IN);


    open(OUT, ">$outhtml_full");
    my $width = 320;
    my $sep_hight = 20 ;

    printf(OUT "<html>\n\n<head>\n\n</head>\n\n<body>\n");
    printf(OUT "<table border=\"1\" cellpadding=\"5\">\n");
    printf(OUT "  <caption>%s</caption>\n", $caption);
    printf(OUT "\n");
    printf(OUT "\n");

    my $nlist = $#list + 1;
    my $ncol = 4;
    for(my $ilist = 0; $ilist < $nlist; $ilist += $ncol){
        my @fig = ();
        my @id = ();
        for(my $icol = 0; $icol < $ncol; $icol++){
            if($ilist + $icol < $nlist){
                $fig[$icol] = $list[$ilist + $icol];
            }
        }
        printf(OUT "  <tr>\n");
        for(my $icol = 0; $icol < $ncol; $icol++){
            if($ilist + $icol < $nlist){
                printf(OUT "    <th>%s</th>\n", $fig[$icol]);
            }
        }
        printf(OUT "  </tr>\n");
        printf(OUT "  <tr>\n");
        for(my $icol = 0; $icol < $ncol; $icol++){
            if($ilist + $icol < $nlist){
                printf(OUT "    <td><a href=\"./%s\"><img src=\"./%s\"", $fig[$icol], $fig[$icol]);
                printf(OUT " alt=\"%s\" width=\"%d\"> </a> </td>\n", $fig[$icol], $width);
            }
        }
        printf(OUT "  </tr>\n");
        printf(OUT "  <tr height=\"%d\"> </tr>\n", $sep_hight);
        printf(OUT "\n");
        printf(OUT "\n");
    }
    printf(OUT "  </tr>\n</table>\n</body>\n</html>");
    close(OUT);

}
