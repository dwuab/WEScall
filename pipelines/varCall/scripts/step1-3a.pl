#!/usr/bin/perl -w

use strict;
use Time::HiRes qw(usleep nanosleep);
use File::Path;
use File::Basename;
use FindBin;
use lib $FindBin::Bin;
use hyunlib qw(%hszchrs initRef forkExecWait);
use gcconfig;

&initRef($ref);
my($type,$chr,$sitePath,$cpu_step1,$cpu_step2,$cpu_step3) = @ARGV;

my $scriptPath = "$FindBin::Bin/../scripts";

my $result = system("perl $scriptPath/step1-detect-and-merge-variants.pl  $chr");
if($result!=0) {die("Error: perl $scriptPath/step1-detect-and-merge-variants.pl  $chr")};

$result = system("make -f out/aux/chr$chr.Makefile -j  $cpu_step1");
if($result!=0) {die("Error: make -f out/aux/chr$chr.Makefile -j  $cpu_step1")};

if($type eq "WES"){
   # Add the sites existing in 1KG3 reference panel but absent in the discover step.  
   $result = system("perl $scriptPath/add1KG3Sites.pl out/aux/union $chr  $scriptPath $sitePath");
   if($result!=0) {die("Error: perl $scriptPath/add1KG3Sites.pl out/aux/union $chr  $scriptPath $sitePath")};
}

$result = system("perl $scriptPath/step2-joint-genotyping.pl  $chr");
if($result!=0) {die("Error: perl $scriptPath/step2-joint-genotyping.pl  $chr")};

$result = system("make -f out/paste/chr$chr.Makefile -j  $cpu_step2");
if($result!=0) {die("Error: make -f out/paste/chr$chr.Makefile -j  $cpu_step2")};

$result = system("perl $scriptPath/step3a-compute-milk-score.pl $chr");
if($result!=0) {die("Error: perl $scriptPath/step3a-compute-milk-score.pl $chr")};

$result = system("make -f out/aux/milk/milk.$chr.$chr.Makefile -j $cpu_step3");
if($result!=0) {die("Error: make -f out/aux/milk/milk.$chr.$chr.Makefile -j $cpu_step3")};

$result = system("touch out/aux/milk/milk.$chr.$chr.Milk.OK");
if($result!=0) {die("Error: touch out/aux/milk/milk.$chr.$chr.Milk.OK")};
