#!/usr/bin/perl -w

use strict;
use Time::HiRes qw(usleep nanosleep);
use File::Path;
use File::Basename;
use FindBin;
use lib $FindBin::Bin;

my ($mypath, $chr, $TopMed, $refPanel)=@ARGV;

opendir(TMP,$mypath) ;
my $match=$chr."_";
my @file = readdir(TMP);
for(my $i=0; $i<@file; $i++){
	if($file[$i]=~/sites.bcf/ && $file[$i]!~/csi/ && $file[$i]!~/OK/){
		if($file[$i]=~/^$match/ && $file[$i]!~/1KG3/){
			my @t = split(/_|\./,$file[$i]);

			# check if bcftools is located at the expected path
			if (!(-e "$TopMed/../bcftools/bcftools")) {die "$TopMed/../bcftools/bcftools not found!"};

			my $cmd="$TopMed/../bcftools/bcftools  concat $mypath/$file[$i]   $refPanel/ALL.chr$chr.phase3.20130502.SNP.biallelic.MAF0.01.sites.vcf.gz   -r $chr:$t[1]-$t[2] -o $mypath/$file[$i].addSite -O b  -a  -d all";	

			my $result = system($cmd);
			if($result!=0) {die("bcftools concat error!")};
	        system("mv $mypath/$file[$i] $mypath/$t[0]_$t[1]_$t[2].No1KG3.sites.bcf");
	        system("mv $mypath/$file[$i].csi  $mypath/$t[0]_$t[1]_$t[2].No1KG3.sites.bcf.csi");
	        system("mv $mypath/$file[$i].addSite $mypath/$file[$i]");
	        system("$TopMed/../bcftools/bcftools index $mypath/$file[$i]");
	        #Record the sites discovered based on sequencing datasets
	        system("$TopMed/../bcftools/bcftools view  $mypath/$t[0]_$t[1]_$t[2].No1KG3.sites.bcf  | grep -v \"#\" | less -S | cut -f 1,2,3,4,5 >$mypath/$t[0]_$t[1]_$t[2].No1KG3.sites.txt");
	    }
	}
}
