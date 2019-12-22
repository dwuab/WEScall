package gcconfig;

use base qw/Exporter/;
use Cwd qw(realpath);
use File::Basename qw(dirname);
use POSIX qw(pow sqrt);
use FindBin;
use YAML::XS 'LoadFile';

## Variables and methods shared across the package
our @EXPORT = qw($regionBatchSize $REF_PATH $time_latency_job  $ref $md5 $bgzip $samtools $bcftools $bamUtil $tabix $index $pedf $out $vt $discoverUnit $genotypeUnit $batchtype $batchopts_step1 $batchopts_step2  $batchopts_step3 
 $omnivcf $hapmapvcf $dbsnp $invNorm $svmlearn $svmclassify $vcfsummary $vcfsummary2 $sampleBatchSize $discover_thread_perBatch $jointcall_thread_perBatch $milk_thread_perBatch);

############################################################
### MODIFY THESE VARIABLES TO YOUR COMPUTING ENVIRONMENT
our $index = "./data/samples.index";
our $pedf = "./data/samples.ped";
our $out = "out";

### clster settings are located at an external yaml configure file
my $config = LoadFile("$FindBin::Bin/../../../cfg/varCall.cfg.yaml") or die __FILE__ . "configure file ../../cfg/varCall.cfg.yaml not found!\n";

our $sampleBatchSize = $config->{"sampleBatchSize"}; 
our $regionBatchSize = $config->{"regionBatchSize"}; 
our $discover_thread_perBatch = $config->{"discover_thread_perBatch"};    
our $jointcall_thread_perBatch = $config->{"jointcall_thread_perBatch"};
our $milk_thread_perBatch = $config->{"milk_thread_perBatch"};
our $discoverUnit = $config->{"discoverUnit"};
our $genotypeUnit = $config->{"genotypeUnit"};
our $time_latency_job = $config->{"time_latency_job"};

our $batchtype = $config->{"batchtype"};
our $batchopts_step1 = $config->{"batchopts_step1"};
our $batchopts_step2 = $config->{"batchopts_step2"};
our $batchopts_step3 = $config->{"batchopts_step3"};

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
