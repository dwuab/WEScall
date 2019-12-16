#!/usr/bin/perl -w

use strict;
use Time::HiRes qw(usleep nanosleep);
use File::Path;
use File::Basename;
use FindBin;
use lib $FindBin::Bin;
use hyunlib qw(%hszchrs initRef);
use gcconfig;

&initRef($ref);

my @chrs = @ARGV;

if ( $#chrs < 0 ) {
    die "Usage: [command] [list of chromosomes separated by space]\n";
}

my $milkDir = "$out/aux/milk";
my $pasteDir = "$out/paste";
my $unit = $genotypeUnit;

unless ( -e "$milkDir" ) {
    mkpath("$milkDir") || die "Cannot create directory $pasteDir\n";
    
}

unless (-e "$milkDir/jobfiles"){
    mkpath("$milkDir/jobfiles") || die "Cannot create directory $pasteDir $milkDir/jobfiles\n";
}


my $outf = "$milkDir/milk.$chrs[0].$chrs[$#chrs]";

open(MAK,">$outf.Makefile") || die "Cannot open file $outf.Makefile for writing\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all:";
foreach my $chr (@chrs) {
    my $szchr = $hszchrs{$chr}->[3];    
    my $catvcf = "$milkDir/$chr\_1\_$szchr\_milk";
    print MAK " $catvcf.sites.vcf.gz.tbi"
}
print MAK "\n\n";

foreach my $chr (@chrs) {
    unless ( -e "$milkDir/$chr" ) {
    mkdir("$milkDir/$chr") || die "Cannot create directory $pasteDir\n";
    }
    
    my $szchr = $hszchrs{$chr}->[3];
    my @cmds = ();
    my @tgts = ();
    my @milkVcfs = ();
    my @milkVcf_Sites = ();

    for(my $i=1; $i < $szchr; $i += $unit) {
    my $beg = $i;
    my $end = ( $i + $unit > $szchr ) ? $szchr : $i + $unit - 1;
    my $pasteBcf = "$pasteDir/$chr/$chr\_$beg\_$end\_paste.bcf";
    my $milkVcf = "$milkDir/$chr/$chr\_$beg\_$end\_milk.vcf.gz";
    my $milkVcf_Site = "$milkDir/$chr/$chr\_$beg\_$end\_milk.sites.vcf.gz";

    my $cmd = "set -o pipefail; $REF_PATH $vt milk_filter -f $pedf -b $pasteBcf -o $milkVcf && $tabix -pvcf $milkVcf && zcat $milkVcf |  cut -f 1-8 | $bgzip -c > $milkVcf_Site && $tabix -p vcf  $milkVcf_Site";
    push(@cmds,$cmd);
    push(@milkVcfs,$milkVcf);
    push(@tgts,"$milkVcf.tbi");
    push(@milkVcf_Sites, $milkVcf_Site);
    }

    my $catvcf = "$milkDir/$chr\_1\_$szchr\_milk";
    my $cmd = getClusterCmd("make -f $milkDir/jobfiles/milk.$chr.Makefile -j $milk_thread_perBatch ");
    print MAK "$catvcf.sites.vcf.gz.tbi: $milkDir/jobfiles/milk.$chr.OK\n";
    print MAK "\t($tabix -h $milkVcf_Sites[0] NA:0; zcat @milkVcf_Sites | grep -v ^#;) | $bgzip -c > $catvcf.sites.vcf.gz\n";
    print MAK "\t$tabix -pvcf $catvcf.sites.vcf.gz\n\n";
    print MAK "\n$milkDir/jobfiles/milk.$chr.OK:\n";
    print MAK "$cmd\n";

    
   # for(my $j=0; $j < @tgts; ++$j) {
   # print MAK "$tgts[$j]:\n";
   # print MAK "\t$cmds[$j]\n\n";
   # } 


    my $cnt = -1;
    open O, ">$milkDir/jobfiles/milk.$chr.Makefile";
    print O ".DELETE_ON_ERROR:\n\n";
    print O "all: @tgts\n\n";
    print O "\ttouch $milkDir/jobfiles/milk.$chr.OK\n";
    for(my $j=0; $j < @tgts; ++$j) {
        print O "$tgts[$j]:\n";
        $cnt ++;
        #my $MyCMD=getClusterCmd($cmds[$j]);
        #print MAK "$MyCMD\n\n";
        print O "\t$cmds[$j]\n\n";
    }  
    close O;
}
close MAK;

print "Run make -f $outf.Makefile -j [numjobs] to complete this step\n";


sub getClusterCmd {
    my ($cmd, $cmdKey) = @_;

    my $logOption = "";
    my $scriptPath=$FindBin::Bin;
    #if(defined ($cmdKey) && $cmdKey)
    #{
        $logOption = "-log $outf.Makefile.cluster ";
    #}

    $cmd =~ s/'/"/g;            # Avoid issues with single quotes in command
    my $runcluster="$scriptPath/runcluster.pl";
    my $outputDir=$out;
    my $newcmd = $runcluster. "  -bashdir $outputDir/paste/jobfiles  ${logOption}";
    my $BatchOpt="$batchopts_step3  -e  $outputDir/log/   -o  $outputDir/log/";
    if($BatchOpt)
    {
        $newcmd .= "-opts '".$BatchOpt."' ";
    }
    $newcmd .= "$batchtype    '$cmd'";

    $newcmd="\t$newcmd";
    return $newcmd;
}
