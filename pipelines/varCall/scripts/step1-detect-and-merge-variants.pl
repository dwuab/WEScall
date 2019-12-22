#!/usr/bin/perl -w

use warnings;
use Time::HiRes qw(usleep nanosleep);
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use hyunlib qw(%hszchrs initRef);
use lib "./scripts";
use gcconfig;

die "Usage: [command] [list of space separated chromosomes]\n" if ( $#ARGV < 0 );

my $batch_size= $sampleBatchSize;
my $outputDir = $out;
my $sleep = 0;
my $sampleFile = $index;
my $intervalWidth = $discoverUnit; #20000000;
my $chromosomes = join(",",@ARGV);
my $makeFile = "$outputDir/aux/chr$chromosomes.Makefile";
my $refGenomeFASTAFile = $ref;
my $generateIntermediateFiles = 0;
my $useClipOverlap = 0;
my $discoverOptions = "-z -q 20";
my $mergeCandidateVariantsOptions = "";
my $partition = "cluster";
my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my @cmd;
my $inputVCFFile;
my $inputVCFFileList;
my $outputVCFFile;
my $processByGenome = ($intervalWidth==0);
mkpath($outputDir);
my $logDir = "$outputDir/log";
mkpath($logDir);
my $auxDir = "$outputDir/aux";
mkpath($auxDir);
my $individualDir = "$outputDir/aux/individual";
mkpath($individualDir);
my $unionDir = "$outputDir/aux/union";
mkpath($unionDir);
my $evaluationDir = "$outputDir/aux/evaluation";
mkpath($evaluationDir);
my $finalDir = "$outputDir/aux/sites";
mkpath($finalDir);
my $slurmScriptsDir = "$outputDir/slurm_scripts";
mkpath($slurmScriptsDir);
mkpath("$individualDir/../jobfiles");
my $slurmScriptNo = 0;
my $logFile = "$logDir/run.log";

###########################
#Read samples and BAM paths
###########################
my %BAMFILE = ();
my @SAMPLE = ();
print "$sampleFile\n";

open(INDEX,"$sampleFile") || die "Cannot open $sampleFile\n";
while (<INDEX>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($sampleID, $bamPath) = split('\t|\s+', $_);
        $BAMFILE{$sampleID} = $bamPath;
        print "$bamPath  $sampleID\n";
        push(@SAMPLE, $sampleID);
        mkpath("$individualDir/$sampleID");
    }
}
close(INDEX);
print STDOUT "@SAMPLE\n";


###################
#Generate intervals
###################
my %intervalNamesByChrom = ();
my @intervalNames = ();
my @intervals = ();
my @CHROM = ();

my $refGenomeFASTAIndexFile = "$refGenomeFASTAFile.fai";
open(SQ," $refGenomeFASTAIndexFile ") || die "Cannot open  $refGenomeFASTAIndexFile \n";
my %CHROM = ();
map {$CHROM{$_}=1} split(",", $chromosomes);
while (<SQ>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($chrom, $len, $rest) = split('\t', $_);

        next if (!exists($CHROM{$chrom}));

        print STDOUT basename(__FILE__) . ": processing $chrom with length: $len bp\n";

        push(@CHROM, $chrom);

        $intervalNamesByChrom{$chrom} = ();
        my $count = 0;
        for my $i (0 .. floor($len/$intervalWidth))
        {
            my $interval = "";
            my $intervalName = "";
            if ($i<floor($len/$intervalWidth))
            {
                $interval = $chrom . ":" . ($intervalWidth*$i+1) . "-" . ($intervalWidth*($i+1));
                $intervalName = $chrom . "_" . ($intervalWidth*$i+1) . "_" . ($intervalWidth*($i+1));
            }
            elsif ($i*$intervalWidth!=$len)
            {
                $interval = $chrom . ":" . ($intervalWidth*$i+1) . "-" . $len;
                $intervalName = $chrom . "_" . ($intervalWidth*$i+1) . "_" . $len;
            }
            else
            {
                last;
            }

            $count++;

            push(@intervals, $interval);
            push(@intervalNames, $intervalName);
            push(@{$intervalNamesByChrom{$chrom}}, $intervalName);
        }

        print STDOUT basename(__FILE__) . "added $count intervals\n";
    }
}
close(SQ);

#############
#1. Discovery
#############

#log start time
$tgt = "$logDir/start.discovery.OK";
$dep = "";
@cmd = "date | awk '{print \"vt calling pipeline\\n\\nstart discovery: \"\$\$0}' > $logFile";
makeLocalStep($tgt, $dep, @cmd);

my @intervalSampleDiscoveryVCFFilesOK = ();

#generate slurm job array script
my $noJobs = scalar(@intervals)*scalar(@SAMPLE);
my $no = 0;
mkpath("$slurmScriptsDir/job_array_output");
my $slurmJobArrayScript = "$slurmScriptsDir/job_array.sh";
open(SCRIPT, ">$slurmJobArrayScript");
print SCRIPT <<SCRIPT;
\#!/bin/bash
\#SBATCH --partition=main,nomosix
\#SBATCH --error=$slurmScriptsDir/job_array_output/%a_%N_%j.err
\#SBATCH --output=$slurmScriptsDir/job_array_output/%a_%N_%j.log
\#SBATCH --job-name=vt_discovery
#\#SBATCH --time=6:0:0
#\#SBATCH --mem=2G
\#SBATCH --cpus-per-task=1
\#SBATCH --array=1-$noJobs
declare -a commands
SCRIPT

#mine variants from aligned reads
my @batch_sample = ();
for my $i (0 .. $#intervals)
{
    my $intervalVCFFilesOK = "";
    my $SampleCnt = 0;
    my @cmd=();
    my @cmd_OK=();
    my $preSampleCnt = 0;
    my $sample_max = @SAMPLE;

    print STDOUT basename(__FILE__) . ": @SAMPLE\n";

    for my $sampleID (@SAMPLE)
    {
        print STDOUT basename(__FILE__) . ": processing" . $sampleID."\n";
        $SampleCnt ++;
        $outputVCFFile = "$individualDir/$sampleID/$intervalNames[$i].sites.bcf";
        #$tgt = "$outputVCFFile.OK";
        $dep = "$logDir/start.discovery.OK";
        my $current_cmd="set -o pipefail; $REF_PATH  $samtools view -h $BAMFILE{$sampleID} $intervals[$i] -u | $bamUtil clipoverlap --poolSize 100000000 --in -.ubam --out -.ubam | $vt discover2 -z -q 20 -b + -r $refGenomeFASTAFile -s $sampleID -i $intervals[$i] -o $outputVCFFile 2> $individualDir/$sampleID/$intervalNames[$i].discover2.log && touch $outputVCFFile.OK";
        push(@cmd,$current_cmd);
        print STDOUT basename(__FILE__) . ": command to be executed: $current_cmd\n";
        push(@cmd_OK,"$outputVCFFile.OK");
        if($SampleCnt % $batch_size == 0 || $SampleCnt == $sample_max ){
            my $batchName = "$intervalNames[$i]"."_S$preSampleCnt"."_$SampleCnt";
            $tgt = "$individualDir/$batchName.OK";
            $intervalVCFFilesOK .= " $tgt";
            my $dep_tmp = join(" ", @cmd_OK);
            open O, ">$individualDir/../jobfiles/$batchName.job.Makefile";
            print O ".DELETE_ON_ERROR:\n\n.PHONY: clean\n";
            print O "\nall: $tgt\n";
            print O "\n$tgt: $dep_tmp\n";
            print O "\ttouch $tgt\n";

            for(my $i=0;$i<@cmd;$i++){
                print O "\n$cmd_OK[$i]:\n";
                print O "\t$cmd[$i]\n";
            }
            close O;
            makeJob($partition, $tgt, $dep, "make -f $individualDir/../jobfiles/$batchName.job.Makefile -j $discover_thread_perBatch  ");
            @cmd=();
            @cmd_OK=();
            $preSampleCnt = $SampleCnt;
        }


        
        #print SCRIPT "commands[" . ++$no . "]= [ ! -e $outputVCFFile.OK ] && $slurmScriptsDir/$slurmScriptNo.sh && touch $outputVCFFile.OK;\n";
        print SCRIPT "commands[" . ++$no . "]= $slurmScriptsDir/$slurmScriptNo.sh && touch $outputVCFFile.OK;\n";

        
    }
    push(@intervalSampleDiscoveryVCFFilesOK, $intervalVCFFilesOK);
}

print SCRIPT "bash -c \"\${commands[\${SLURM_ARRAY_TASK_ID}]}\"\n";
close(SCRIPT);

for my $i (0..$#intervals)
{
    $inputVCFFileList = "$unionDir/$intervalNames[$i]_vcf_file.list";
    open(OUT, ">$inputVCFFileList");
    for my $sample (@SAMPLE) {print OUT "$individualDir/$sample/$intervalNames[$i].sites.bcf\n";}
    close(OUT);

    #merge variants and annotate variants
    $outputVCFFile = "$unionDir/$intervalNames[$i].sites.bcf";
    $tgt = "$outputVCFFile.OK";
    $dep = $intervalSampleDiscoveryVCFFilesOK[$i];
    @cmd = ("$vt merge_candidate_variants2 $mergeCandidateVariantsOptions -L $inputVCFFileList -o + 2> $unionDir/$intervalNames[$i].merge_candidate_variants2.log | " .
            "$vt annotate_indels -r $refGenomeFASTAFile + -o + 2> $unionDir/$intervalNames[$i].annotate_indels.log | " .
            "$vt consolidate_variants + -o $outputVCFFile 2> $unionDir/$intervalNames[$i].consolidate_variants.log");
    makeJob($partition, $tgt, $dep, @cmd);

    $inputVCFFile = "$unionDir/$intervalNames[$i].sites.bcf";
    $tgt = "$inputVCFFile.csi.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt index $inputVCFFile");
    makeJob($partition, $tgt, $dep, @cmd);
}



#log end time
$tgt = "$logDir/end.discovery.OK";
$dep = join(" ", map {"$unionDir/$_.sites.bcf.OK"} @intervalNames);
@cmd = ("date | awk '{print \"end discovery: \"\$\$0}' >> $logFile");
makeJob("local", $tgt, $dep, @cmd);


##############
#2. Genotyping
##############

#my @intervalSampleGenotypingVCFFilesOK = ();

#mine variants from aligned reads
#for my $i (0 .. $#intervals)
#{
#    my $intervalVCFFilesOK = "";
#    for my $sampleID (@SAMPLE)
#    {
#        $inputVCFFile = "$unionDir/$intervalNames[$i].sites.bcf";
#        $outputVCFFile = "$individualDir/$sampleID/$intervalNames[$i].genotypes.bcf";
#        $tgt = "$outputVCFFile.OK";
#        $dep = "$inputVCFFile.OK";
#  #        @cmd = ("REF_PATH=/dept/csg/topmed/working/mktrost/gotcloud.ref/md5/%2s/%2s/%s; $vt genotype2 $inputVCFFile -b $BAMFILE{$sampleID} -r $refGenomeFASTAFile -s $sampleID -i $intervals[$i] -o $outputVCFFile 2> $individualDir/$sampleID/$intervalNames[$i].genotype2.log");
#        @cmd = ("$vt genotype2 $inputVCFFile -b $BAMFILE{$sampleID} -r $refGenomeFASTAFile -s $sampleID -i $intervals[$i] -o $outputVCFFile 2> $individualDir/$sampleID/$intervalNames[$i].genotype2.log");
#        makeJob($partition, $tgt, $dep, @cmd);
#
#        $inputVCFFile = "$individualDir/$sampleID/$intervalNames[$i].genotypes.bcf";
#        $tgt = "$inputVCFFile.csi.OK";
#        $dep = "$inputVCFFile.OK";
#        @cmd = ("$vt index $inputVCFFile");
#        makeJob($partition, $tgt, $dep, @cmd);
#
#        $intervalVCFFilesOK .= " $outputVCFFile.OK";
#    }
#
#    push(@intervalSampleGenotypingVCFFilesOK, $intervalVCFFilesOK);
#}
#
#my $intervalSampleGenotypingVCFFiles = join(" ", @intervalSampleGenotypingVCFFilesOK);

#log end time
#$tgt = "$logDir/end.genotyping.OK";
#$dep = $intervalSampleGenotypingVCFFiles;
#@cmd = ("date | awk '{print \"end genotyping: \"\$\$0}' >> $logFile");
#makeJob("local", $tgt, $dep, @cmd);

####################
#Write out make file
####################

open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK ".PHONY: clean\n\n";
print MAK "all: @tgts\n\n";

#clean
$tgt = "clean";
$dep = "";
@cmd = ("-rm -rf $finalDir/*.OK $individualDir/*.OK $individualDir/*/*.OK $logDir/*.OK");
makePhonyJob($tgt, $dep, @cmd);
my $cnt = -1;
for(my $i=0; $i < @tgts; ++$i)
{
    
    print MAK "$tgts[$i]: $deps[$i]\n";
    if($tgts[$i]=~/indiv/){
            $cnt ++;
            my $sleepTime = $cnt*$time_latency_job;
            print MAK "\tsleep $sleepTime && $cmds[$i]\n";
    }
    else{
                print MAK "$cmds[$i]\n";
    }
}
close MAK;

print "Run make -f $makeFile -j [numjobs] to complete this step\n";

##########
#functions
##########

#run a job either locally or by slurm
sub makeJob
{
    my ($method, $tgt, $dep, @cmd) = @_;

    if ($method eq "local")
    {
        makeLocalStep($tgt, $dep, @cmd);
    }
    else
    {
        #makeSlurm($partition, $tgt, $dep, @cmd);
        makeClusterStep($tgt, $dep, @cmd);
    }
}

#run slurm jobs
sub makeSlurm
{
    my ($partition, $tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        #contains pipe
        if ($c=~/\|/)
        {
            ++$slurmScriptNo;
            my $slurmScriptFile = "$slurmScriptsDir/$slurmScriptNo.sh";
            open(IN, ">$slurmScriptFile");
            print IN "#!/bin/bash\n";
            print IN "set -o pipefail; $c";
            close(IN);
            chmod(0755, $slurmScriptFile);

            $cmd .= "\techo '" . $c . "'\n";
            #$cmd .= "\tsrun -p $partition $slurmScriptFile\n";
            $cmd .= "\tbash $slurmScriptFile\n";	    
        }
        else
        {
            #$cmd .= "\tsrun -p $partition " . $c . "\n";
            $cmd .= "\tbash " . $c . "\n";	    
        }
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

#run a local job
sub makeLocalStep
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\tset -o pipefail; " . $c . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

#run a local phony job
sub makePhonyJob
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . $c . "\n";
    }
    push(@cmds, $cmd);
}



#run a cluster job
sub makeClusterStep
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        #$cmd .= "set -o pipefail; " . $c . ";";
        $cmd .= "set -o pipefail; " . $c;                           # Note I change this cmd...
    }
    $cmd .= "&& touch $tgt";
    $cmd=getMosixCmd($cmd);
    push(@cmds, $cmd);
}



sub getMosixCmd {
    my ($cmd, $cmdKey) = @_;

    my $logOption = "";
    my $makef_OUT_DIR="$outputDir";
    my $scriptPath=$FindBin::Bin;
    #if(defined ($cmdKey) && $cmdKey)
    #{
    $logOption = "-log $outputDir/aux/chr$chromosomes.Makefile.cluster ";
    #}

    $cmd =~ s/'/"/g;            # Avoid issues with single quotes in command
    my $runcluster="$scriptPath/runcluster.pl";
    my $newcmd = $runcluster. "  -bashdir $outputDir/aux/jobfiles  ${logOption}";
    my $BatchOpt="$batchopts_step1  -e  $outputDir/log/   -o  $outputDir/log/";     
    if($BatchOpt)
    {
        $newcmd .= "-opts '".$BatchOpt."' ";
    }
    $newcmd .= "$batchtype    '$cmd'";
    $newcmd="\t$newcmd";
    return $newcmd;
}
