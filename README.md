# WEScall

#### Authors: Jinzhuang Dou, [Chaolong Wang](http://chaolongwang.github.io)

#### License: [GNU General Public License v3.0 (GPLv3)](https://www.gnu.org/licenses/gpl-3.0.en.html)
---

## 1. Description

WES genotype calling (WEScall) pipeline focuses on variant calling of whole exome sequencing (WES) data. It can call genotypes accurately in both target and off-target regions by making full utilization of off-target reads. 
WEScall could: 
* work on the HPC cluster, supporting both PBS and SGE environments.
* handle cluster specifics internally and users don't worry about the scheduler usage.
* run automatically in parallel to make optimal use of resources where possible. 

## 2. Citation for our pipeline 

Details of the variant calling pipeline can be found in our paper:  
* Dou J, Wu D, ... , Wang C. Joint analysis of target and off-target data improves whole-exome sequencing studies (in preparation)

## 3. Dependencies
* python (version >= 3.5)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (version >= 3.9.1)
* java (version >= 1.8.0)
* perl (version >= v5.10)

## 4. Download and install

You can download our pipeline by typing the following command in a terminal.

`git clone https://github.com/jinzhuangdou/WEScall.git` 

This command will create a folder named "WEScall". 

## 5. Getting started 

**Shorthands**
* **PL_DIR**: the directory where the pipeline is installed. 
* **WK_DIR**: the directory where you store the outputs.

### 5.1. Variant calling 
Before running WEScall, you should first prepare for the sequence index and pedigree files. These two files have the same format required [TopMed](https://github.com/statgen/topmed_freeze3_calling) pipeline.
Each line in the index file contains:
```
  [sampleID] [Full Path to BAM/CRAM file] [Contamination rate -- set zero if unknown].
``` 
Please make sure that all the BAM/CRAM files have been indexed using samtools. 
Each line in the pedigree file includes 
```
  [family ID] [sample ID] [father ID] [mother ID] [gender] [phenotype].
``` 
If no pedigree/ phenotype information available, please set zero in this file. Both the index and pedigree files are tab-delimited. Then you should prepare for the configure file specifying resources required by the pipeline. One example `user.cfg.yaml` in `${PL_DIR}/example/test_WES` is as following: 

```
  chrs: 20-22   
  targetBed:  ${PL_DIR}/WEScall/resources/SeqCap_EZ_Exome_v3_primary.bed
  1KG3_panel: ${PL_DIR}/WEScall/resources/data_v5a_filtered
  geneticMap: ${PL_DIR}/WEScall/resources/geneticMap_GRCh37
``` 
The first line specifies which chromosomes you want to call. Chromosomes that can be called are chromosomes 1 to 22 and X. WEScall coud not call variants from chromosome Y. You can specify multiple chromosomes one by one, delimited by dash (for example 20-22). 

The second line specifies the target region bed file. It lists the targeted exonic regions with start and stop chromosome locations in GRCh37/hg19. Note, WEScall can also support the analysis of samples from whole genome sequencing (WGS). In this case, the target region bed file is not necessary. 

The third line specifies the location of the 1KG3 reference panel used for genotype imputation. 

The fourth line specifies the location of the genetic map files used for genotype imputation. 

Now we can generate the master job file using the following command

```
  python  ${PL_DIR}/varCall/runTopMed.py -c  user.cfg.yaml   -i samples.index -p samples.ped  -t seq_type WES -n 
``` 
The `seq_type`option has to be either WGS or WES. WEScall can also analyze samples from whole genome sequencing (WGS). After running this command, there will generate the folder `${WK_DIR}/varCall` storing the execute script `${WK_DIR}/varCall/run.sh` and configure file `${WK_DIR}/varCall/cluster.yaml`. Users can modify these files before running the pipeline if necessary. 

You can submit the variant calling master job using
```
  cd varCall & qsub run.sh >> ./logs/submission.log  
``` 
If any job is killed prematurely, you can resume the master job by using the command again. You can check `${WK_DIR}/varCall/logs/snakemake.log` for progress or diagnose premature terminations of jobs. Once snakemake.log reports all jobs are done, you can proceed to the next step.

### 5.2. Divide the large genome into short chromosomal segments 

A good way to avoid the huge memory usage in downstream imputation procedure is to split the genome into smaller regions. Do this by 
running the following command:
```
  python ${PL_DIR}/phasing/splitGenome.py  -t seq_type -c user.cfg.yaml -n  
``` 
After running this command, there will generate the folder `${WK_DIR}/phasing` storing the script `${WK_DIR}/phasing/run.sh` and users can submit the genome splitting master job using
```
  cd phasing
  qsub run.sh >> ./logs/submission.log
```

### 5.3. Genotype refinement through phasing

This step performs genotype refinement through phasing by leveraging linkage disequilibrium (LD) information from study samples of external reference panel. Run the following command to generate the job file: 
```
  python  ${PL_DIR}/phasing/phasing.py  -t seq_type -c user.cfg.yaml -n
```
Before running the phasing procedure, you can modify parameters in the configure file `${WK_DIR}/varCall/cluster.yaml` if necessary. If you want to change the number of segments being phased at one time, please modify the option `N_ARG="--jobs 200"` in line 81 of file `${WK_DIR}/phasing/run.sh`. 

Finally, you can submit the variant phasing master job using
```
  cd phasing & qsub run.sh >> ./logs/submission.log
```
When all above jobs are finished, the genotyping results are stored in `${WK_DIR}/phasing`. For example, users can see genotypes from chromosome 1 in `${WK_DIR}/phasing/1/1.Final.vcf.gz` 

## 6. Frequently used settings and operations

### 6.1. PBS specific settings
If you want to run this pipeline on the cluster of PBS and require to add the project IDs and queue names. You can do the following modifications.
```
  ${PL_DIR}/pipelines/lib/run.template.PBS.sh    Line:44-45
  44  #PBS -P 13000026
  45  #PBS -q production

  ${PL_DIR}/topMed/scripts/gcconfig.pm Line:26,27,28
  26  our $batchopts_step1 = "   -l  select=1:ncpus=1:mem=10G  -l walltime=12:00:00 -P 13000026 -q production";
  27  our $batchopts_step2 = "   -l  select=1:ncpus=1:mem=50G -l walltime=24:00:00 -P 13000026 -q production";
  28  our $batchopts_step3 = "   -l  select=1:ncpus=1:mem=50G -l walltime=24:00:00 -P 13000026 -q production";
```
### 6.2. Memory settings of variant calling
The joint calling step may take huge memory when the sample size is very large (>1,000). This may kill the program. You can address this issue by either setting the maximum memory usage or splitting genome into smaller regions (default 1Mb).

Increase the maximum memory usage:
```
${PL_DIR/topMed/scripts/gcconfig.pm Line:27
27  our $batchopts_step2 = "   -l  select=1:ncpus=1:mem=50G -l walltime=24:00:00 -P 13000026 -q production";
```
Split the whole genome into smaller regions (default 1Mb): 
```
${PL_DIR}/topMed/scripts/gcconfig.pm Line:19
19  our $genotypeUnit = 1000000;
```
## 7. Questions
For further questions, pleast contact Jinzhuang Dou <jinzhuangdou198706@gmail.com>.
