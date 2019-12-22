#!/usr/bin/env perl

use strict;
use warnings;
use English;
use Getopt::Long;
use FindBin qw($Bin); 
use lib "$FindBin::Bin";


my $cmdPath = $0;
my $cmd = "";
my $invcf = "";
my $out = "";
my $rmLst = "";
my $keepLst = "";
my $type = "";
my $cmdDir = $cmdPath; $cmdDir =~ s/\/[^\/]+$//; $cmdDir =~ s/\/[^\/]+$//;
my $bgzip = "$FindBin::Bin/htslib/bgzip";

my $result = GetOptions(
    "invcf=s",\$invcf,
	"removeLst=s",\$rmLst,
	"keepLst=s", \$keepLst,
	"type=s", \$type,
	"out=s",\$out,
	);

my $usage = <<END;
--------------------------------------------------------------------------------

apply_SVM_filter.pl: Remove sites that do not pass the SVM filter

Usage: perl apply_SVM_filter.pl --invcf [input VCF file] --out [output Dir] \
	Options:
	--invcf      Input VCF file with rich set of INFO fields
	--removeLst  Input site file with the SVM INFO field
	--keepLst Input site file of 1KG3 reference panel
	--type 	  Data type
	--out        The output prefix
	
--------------------------------------------------------------------------------
END

unless (($result) && (($invcf) || ($out) || ($rmLst) || ($keepLst))) {
	die "Error in parsing options\n$usage\n";
}

my(%target, %target_fail, %keepLst);
########################################################
if($type eq "WGS"){

	if ($rmLst =~ m/\.gz$/){open(IN,"gunzip -c $rmLst |");}
	else{open(IN,$rmLst) || die "Cannot open $rmLst\n";}
	while(<IN>){
		chomp;
		# if($_=~/#/ || $_=~/overlap/){;}  # remove the overlap flag 
		if($_=~/#/){;}  
		else{
			my @f=split(/\s+|\t/,$_);
			my $id = "$f[0]-$f[1]-$f[3]-$f[4]";
			$target{$id}++;
			if($f[6]=~/SVM/ || $f[6]=~/EXHET/){
				$target_fail{$id}++;
			}
	    }
	}
	close(IN);

	open O, "| $bgzip -c >  $out";
	if ($invcf =~ m/\.gz$/){open(IN,"gunzip -c $invcf |");}
	else{open(IN,$invcf) || die "Cannot open $invcf\n";}
	while(<IN>){
		chomp;
		my @f=split(/\s+|\t/,$_);
		if($_=~/#/){print O "$_\n";}
		elsif(exists $target_fail{"$f[0]-$f[1]-$f[3]-$f[4]"}){;}
		else{
			print O "$_\n";
		}
	}
	close(IN);
	close (O);
}
elsif($type eq "WES"){
	my $targetTotal=0;
	my $targetFail=0;
	my $offtarget_1KG3=0;
	if ($rmLst =~ m/\.gz$/){
		if (! -e $rmLst) {
			die "$rmLst does not exists!\n";
		}
		open(IN,"gunzip -c $rmLst |") || die "Cannot open $rmLst\n";}
	else {
		open(IN,$rmLst) || die "Cannot open $rmLst\n";}
	while(<IN>){
		chomp;
		# if($_=~/#/ || $_=~/overlap/){;}  # remove the overlap flag 
		if($_=~/#/){;}  
		else{
			my @f=split(/\s+|\t/,$_);
			my $id = "$f[0]-$f[1]-$f[3]-$f[4]";
			$target{$id}++;
			$targetTotal++;
			if($f[6]=~/SVM/){
				$target_fail{$id}++;
				$targetFail++;
			}
                }
	}
	close(IN);

	print "$rmLst  total markers $targetTotal ... svm filter markers $targetFail\n";
	if ($keepLst =~ m/\.gz$/){
		if (! -e $keepLst) {
			die "$keepLst does not exists!\n";
		}
		open(IN,"gunzip -c $keepLst |") || die "Cannot open $keepLst\n";
	}
	else {
		open(IN,$keepLst) || die "Cannot open $keepLst\n";}
	while(<IN>){
		chomp;
		if($_=~/#/){;}  
		else{
			my @f=split(/\s+|\t/,$_);
			my $id = "$f[0]-$f[1]-$f[3]-$f[4]";
			$keepLst{$id}++;
			$offtarget_1KG3++;
	    }
	}
	close(IN);
	print "$keepLst total markers $offtarget_1KG3\n";


	my $targetInCnt=0;
	my $offtargetInCnt=0;
	open O, "| $bgzip -c >  $out";
	if ($invcf =~ m/\.gz$/){
		if (! -e $invcf) {
			die "$invcf does not exists!\n";
		}
		open(IN,"gunzip -c $invcf |") || die "Cannot open $invcf\n";}
	else {
		open(IN,$invcf) || die "Cannot open $invcf\n";}
	while(<IN>){
		chomp;		
		if($_=~/#/){print O "$_\n";}
		else{
			my @f=split(/\s+|\t/,$_);
			my $id = "$f[0]-$f[1]-$f[3]-$f[4]";
			my $dep_ave = $1 if ( $f[7] =~ /AVGDP=([^;]+)/ );
		     	if($dep_ave>0.1){
				if(exists $target{$id}){
					$targetInCnt++;
					if(exists $keepLst{$id}){print O "$_\n";}
					else{
						if(exists $target_fail{$id}){;}
						else{
							print O "$_\n";
						}
					}
				}
			else{
				### MAF>1% FILTER 
				if(exists $keepLst{$id}){
					$offtargetInCnt++;
					my $af = $1 if ( $f[7] =~ /AF=([^;]+)/ );
					#if($af > 0.01 && $af < 0.99){
						print O "$_\n";
					#}
				}
			}
 		    }
		}
	}
	print "Study samples in target region [$targetInCnt] ... in off-target region [$offtargetInCnt]\n";
	close(IN);
	close (O);

}
else{
	print "Please specify the correct data type : WES|WGS\n";
	exit(-1);
}

