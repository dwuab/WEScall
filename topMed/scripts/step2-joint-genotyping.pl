#!/usr/bin/perl -w

use strict;
use Time::HiRes qw(usleep nanosleep);
use File::Path;
use File::Basename;
use FindBin;
use lib $FindBin::Bin;
use hyunlib qw(%hszchrs initRef);
use lib "./scripts";
use gcconfig;

&initRef($ref);

my @chrs = @ARGV;

if ( $#chrs < 0 ) {
    die "Usage: [command] [list of chromosomes separated by space]\n";
}

my $outDir = "$out/paste";
my $unionDir = "$out/aux/union";

mkpath($outDir) unless ( -e $outDir ); # || die "Cannot! create directory $outDir";
mkpath("$outDir/jobfiles") unless ( -e "$outDir/jobfiles" ); # || die "Cannot! create directory $outDir";

my $outf = "$outDir/chr".join("_",@chrs);
open(MAK,">$outf.Makefile") || die "Cannot open file $outf.Makefile for writing\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all:";
foreach my $chr (@chrs) {
    print MAK " $outDir/$chr.OK";
}
print MAK "\n\n";

my $unit = $genotypeUnit;

foreach my $chr (@chrs) {
    unless ( -e "$outDir/$chr" ) {
	mkdir("$outDir/$chr") || die "Cannot create directory $outDir/$chr\n";
    }

    my $szchr = $hszchrs{$chr}->[3];
    my @cmds = ();
    my @tgts = ();
    my @deps = ();    
    my @outvcfs = ();

    for(my $iD=1; $iD < $szchr; $iD += $discoverUnit) {
	my $begD = $iD;
	my $endD = ( $iD + $discoverUnit > $szchr ) ? $szchr : ($iD + $discoverUnit - 1);
	my $sitevcf = "$unionDir/$chr\_$begD\_$endD.sites.bcf";
	


    my $region = 0; 
    my $preBed = 0;
    my $nowBed = 0;
    my $max_region = 0;
    my @region_cmds = ();
    my @region_tgts = ();
    my @region_deps = ();

    for(my $i=$begD; $i < $endD; $i += $unit) { $max_region ++;}

	for(my $i=$begD; $i < $endD; $i += $unit) {
	    my $beg = $i;
	    my $end = ( $i + $unit > $szchr ) ? $szchr : $i + $unit - 1;
	    my $outvcf = "$outDir/$chr/$chr\_$beg\_$end\_paste.bcf"; 

	    #my $cmd = "REF_PATH=$md5 $vt joint_genotype_sequential -r $ref -L $index -i $chr:$beg-$end -o $outvcf $sitevcf 2> $outvcf.log && $vt index $outvcf && touch $outvcf.OK"; 
	    my $cmd = "set -o pipefail; $REF_PATH  $vt joint_genotype_sequential -r $ref -L $index -i $chr:$beg-$end -o $outvcf $sitevcf 2> $outvcf.log && $vt index $outvcf && touch $outvcf.OK"; 
        $region ++;
        push(@region_cmds,$cmd);
        push(@region_tgts,"$outvcf.OK");
        push(@region_deps,"$sitevcf");
        $nowBed = $end;
        if($region % $regionBatchSize==0 || $region == $max_region){
            my $regionFlag = "$chr\_1\_$szchr"."_R$preBed"."_$nowBed";
            open O, ">$outDir/jobfiles/$regionFlag.Makefile";
            print O ".DELETE_ON_ERROR:\n";
            print O "\nall: $outDir/$chr/$regionFlag.OK\n";
            print O "\n$outDir/$chr/$regionFlag.OK:";
            for(my $i=0; $i<@region_tgts; $i++){
                print O " $region_tgts[$i]";
            }
            print O "\n";
            print O "\ttouch $outDir/$chr/$regionFlag.OK\n";
            for(my $i=0; $i<@region_tgts; $i++){
                print O "\n$region_tgts[$i]: $region_deps[$i].OK\n";
                print O "\t$region_cmds[$i]\n";
            }
            print O "\n";
            close O;
            @region_cmds = ();
            @region_tgts = ();
            @region_deps = ();
            $preBed = $nowBed;

            push(@cmds,"make -f $outDir/jobfiles/$regionFlag.Makefile -j $jointcall_thread_perBatch ");
            push(@tgts,"$outDir/$chr/$regionFlag.OK");

            
            #push(@tgts,"$outvcf.OK");
            push(@deps,"$sitevcf");

        }
        push(@outvcfs,$outvcf);

	}
    }
    my $svcf = "$outDir/$chr\_1\_$szchr\_paste.sites.vcf.gz";

    print MAK "$outDir/$chr.OK : @tgts\n";
    print MAK "\t($bcftools view -G $outvcfs[0]";
    for(my $i=1; $i < @outvcfs; ++$i) {
	print MAK "; $bcftools view -G -H $outvcfs[$i]";
    }
    print MAK ") | $bgzip -c > $svcf\n";
    print MAK "\t$tabix -pvcf $svcf\n";    
    print MAK "\ttouch $outDir/$chr.OK\n\n";
    my $cnt = -1;
    for(my $j=0; $j < @tgts; ++$j) {
           $cnt ++; 
           my $sleepTime = $cnt*$time_latency_job;
	       print MAK "$tgts[$j]: $deps[$j].OK\n";
           my $MyCMD=getClusterCmd($cmds[$j]);
	       print MAK "\tsleep $sleepTime && $MyCMD\n\n";
    }
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
    my $newcmd = $runcluster. "   -bashdir $outputDir/paste/jobfiles  ${logOption}";
    my $BatchOpt="$batchopts_step2  -e  $outputDir/log/   -o  $outputDir/log/";
    if($BatchOpt)
    {
        $newcmd .= "-opts '".$BatchOpt."' ";
    }
    $newcmd .= "$batchtype    '$cmd'";
    $newcmd="\t$newcmd";
    return $newcmd;
}
