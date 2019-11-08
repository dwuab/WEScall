package gcconfig;

use base qw/Exporter/;
use Cwd qw(realpath);
use File::Basename qw(dirname);
use POSIX qw(pow sqrt);
use FindBin;

## Variables and methods shared across the package
our @EXPORT = qw($regionBatchSize $REF_PATH $time_latency_job  $ref $md5 $bgzip $samtools $bcftools $bamUtil $tabix $index $pedf $out $vt $discoverUnit $genotypeUnit $batchtype $batchopts_step1 $batchopts_step2  $batchopts_step3 
 $omnivcf $hapmapvcf $dbsnp $invNorm $svmlearn $svmclassify $vcfsummary $vcfsummary2 $sampleBatchSize $discover_thread_perBatch $jointcall_thread_perBatch $milk_thread_perBatch);

############################################################
### MODIFY THESE VARIABLES TO YOUR COMPUTING ENVIRONMENT
our $index = "./data/test.index";
our $pedf = "./data/test.ped";
our $out = "out";
# If you want to reduce the number of jobs in the cluster, please increase the sampleBatchSize and genotypeUnit 
our $sampleBatchSize = 100; 
our $regionBatchSize = 10; 
our $discover_thread_perBatch = 12;    
our $jointcall_thread_perBatch = 12;
our $milk_thread_perBatch = 18;
our $discoverUnit = 200000000000;
our $genotypeUnit = 1000000;
our $time_latency_job = 10;



### Please modify the cluster configuration informaiton for large-scale datasets!
### the following settings are provided for reference
### For jobs in NSCC, please specify the jobIDs 
# if (-d '/home/users/astar/gis/userrig'){
# 	our $batchtype = "pbs";
# 	our $batchopts_step1 = "   -l  select=1:ncpus=$discover_thread_perBatch:mem=24G  -l walltime=24:00:00 -P 13000026 -q production";
# 	our $batchopts_step2 = "   -l  select=1:ncpus=$jointcall_thread_perBatch:mem=32G -l walltime=24:00:00 -P 13000026 -q production";
# 	our $batchopts_step3 = "   -l  select=1:ncpus=$milk_thread_perBatch:mem=24G -l walltime=24:00:00 -P 13000026 -q production";

# }
# elsif(-d '/home/userrig'){
# 	our $batchtype = "sge";
# 	our $batchopts_step1 = "-q medium.q  -pe OpenMP 1  -l mem_free=4G,h_rt=12:00:00 -V -cwd -terse -b y";
# 	our $batchopts_step2 = "-q medium.q  -pe OpenMP 1  -l mem_free=12G,h_rt=12:00:00 -V -cwd -terse -b y";
# 	our $batchopts_step3 = "-q medium.q  -pe OpenMP 1  -l mem_free=12G,h_rt=48:00:00 -V -cwd -terse -b y";
# }
# else{
# 	our $batchtype = "pbs";
# 	our $batchopts_step1 = "-q cu -l mem=4g -V";
# 	our $batchopts_step2 = "-q cu -l mem=12g -V";
# 	our $batchopts_step3 = "-q cu -l mem=12g -V";
# }

our $batchtype = "pbs";
our $batchopts_step1 = "-q cu -l mem=4g -V";
our $batchopts_step2 = "-q cu -l mem=12g -V";
our $batchopts_step3 = "-q cu -l mem=12g -V";

############################################################
### MODIFY THESE VARIABLES TO IF REFERENCE IS LOCATED ELSEWHERE
our $REF_PATH = "";  
our $refDir = "$FindBin::Bin/../gotcloud.ref";
our $md5 = "$refDir/md5/%2s/%s/%s";
our $ref = "$refDir/hs37d5.fa";
our $dbsnp = "$refDir/dbsnp_142.b37.vcf.gz";
our $hapmapvcf = "$refDir/hapmap_3.3.b37.sites.vcf.gz";
our $omnivcf = "$refDir/1000G_omni2.5.b37.sites.PASS.vcf.gz";


############################################################
### MODIFY THESE VARIABLES TO IF EXTERNAL BINARIES ARE USED
our $bgzip = "$FindBin::Bin/../htslib/bgzip";
our $tabix = "$FindBin::Bin/../htslib/tabix";
our $vt = "$FindBin::Bin/../vt/vt";
our $samtools = "$FindBin::Bin/../samtools/samtools";
our $bcftools = "$FindBin::Bin/../bcftools/bcftools";
our $bamUtil = "$FindBin::Bin/../gotcloud/src/bin/bamUtil";
our $svmlearn = "$FindBin::Bin/../gotcloud/src/bin/svm-train";
our $svmclassify = "$FindBin::Bin/../gotcloud/src/bin/svm-predict";
our $invNorm = "$FindBin::Bin/../gotcloud/src/bin/invNorm";
our $vcfsummary = "$FindBin::Bin/vcf-summary";
our $vcfsummary2 = "$FindBin::Bin/vcf-summary-v2";
1;
