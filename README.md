# WEScall

#### Authors: Jinzhuang Dou, Degang Wu, [Chaolong Wang](http://chaolongwang.github.io)

#### License: [GNU General Public License v3.0 (GPLv3)](https://www.gnu.org/licenses/gpl-3.0.en.html)
---

## 1. Description

The WES genotype calling (WEScall) pipeline can perform variant calling of whole exome sequencing (WES) data as well as whole genome sequencing (WGS) data. It can call genotypes accurately in both target and off-target regions by making full utilization of off-target reads. 
WEScall can: 
* runs on HPC cluster, can be configured to support PBS, SGE and other environments.
* handles cluster specifics internally and users don't need to worry about the scheduler usage.
* runs automatically in parallel to make optimal use of resources where possible. 

## 2. Citation for our pipeline 

Details of this pipeline can be found in our paper:  
* Dou J, Wu D, ... , Wang C. Joint analysis of target and off-target data improves whole-exome sequencing studies (in preparation)

## 3. Dependencies
* python (version >= 3.5)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (version >= 3.9.1)
* java (version >= 1.8.0)
* perl (version >= v5.10)
* bcftools (version >= 1.9)

## 4. Download and install

You can download our pipeline by the following command:

`git clone https://github.com/dwuab/WEScall.git` 

This command will create a folder named "WEScall". 

## 5. Getting started 

**Environment variables**
* **PL_DIR**: the directory where the pipeline is located. 
* **WK_DIR**: the directory where you store the outputs.

### 5.1 Setting up for your cluster

As there are no reliable ways to automatically detect the type of job scheduling systems (SGE, PBS, etc.), you have to do some set-up before running the pipeline. You should check `${PL_DIR}/lib/run.template.sh` and `${PL_DIR}/topMed/scripts/gcconfig.pm` and modify them if necessary. The default settings are tested on a Torque (An implementation of PBS) cluster. Earlier versions of this pipeline have been tested on a SGE cluster and a PBS Pro cluster. Template files for SGE and PBS Pro clusters are provided for reference under `${PL_DIR}/lib/`. Settings for SGE and PBS Pro clusters are also provided in the comments in `${PL_DIR}/topMed/scripts/gcconfig.pm`.

### 5.2 Generating 1KG reference panel

**Please run** `create_g1k_ref.sh` to generate the files needed from official 1KGP3 release .vcf files.

### 5.3 Setting up links to resource files

**Please run** `check_resources.sh` to check what resources files you lack and ways to download it.

## 6. Running the pipeline

### 6.1. Variant calling 
Before running WEScall, you should first prepare a file, `samples.index`, containing a list of samples to call. The file have the same format required by [TopMed](https://github.com/statgen/topmed_freeze3_calling) pipeline.
Each line of `samples.index` is of the following format:
```
  [sampleID] [Full Path to BAM/CRAM file] [Contamination rate -- set to zero if unknown].
``` 
**THERE CANNOT BE EMPTY LINES IN `samples.index` FILE!**
The index file has to be **tab-delimited**. 
BAM/CRAM files listed are assumed to be indexed and contain no hard clipped reads, i.e., reads whose CIGAR string contains "H".

You should also prepare a configure file specifying the chromosomes to call, the paths to the resources required by the pipeline and the type of sequencing data (WES or WGS). One example `user.cfg.yaml` in `${PL_DIR}/example/test_WES` is as following: 

```
  chrs: 1,2,10,20,X  
  targetBed:  ${PL_DIR}/WEScall/resources/SeqCap_EZ_Exome_v3_primary.bed
  1KG3_panel: ${PL_DIR}/WEScall/resources/data_v5a_filtered
  geneticMap: ${PL_DIR}/WEScall/resources/geneticMap_GRCh37
  seqType:    WES
``` 
The first line specifies which chromosomes you want to call. Chromosomes that can be called are chromosomes 1 to 22 and X. WEScall can not call variants from chromosome Y. You can specify multiple chromosomes one by one, delimited by comma (for example 20,22). 

The second line specifies the target region bed file. It lists the targeted exonic regions with start and stop chromosome locations in GRCh37/hg19. Note, WEScall can also support the analysis of WGS samples, in which case the target region bed file is not necessary. 

The third line specifies the location of the 1KG3 reference panel used for genotype phasing. Please see the citation for the filtering criteria.

The fourth line specifies the location of the genetic map files used for genotype phasing. 

The fifth line specifies the type of the sequence data. The allowable values are WES and WGS.

Now we can generate the master job file using the following command

```
  cd ${WK_DIR} && python ${PL_DIR}/WEScall.py varCall -c user.cfg.yaml -s samples.index
``` 
After running this command, the folder `${WK_DIR}/varCall` will be generated, storing the execute script `${WK_DIR}/varCall/run.sh` and configure file `${WK_DIR}/varCall/cluster.yaml`. Users can modify these files before running the pipeline if necessary. 

You can submit the variant calling master job using
```
  cd varCall & qsub run.sh >> ./logs/submission.log  
``` 
If any job is killed prematurely, you can resume the master job by using the command again. You can check `${WK_DIR}/varCall/logs/WEScall_varCall.master.log` for progress or diagnose premature terminations of jobs. Once the log reports all jobs are done (message such as "4 of 4 steps (100%) done"), you can proceed to the next step.

### 6.2. Divide the large genome into short chromosomal segments 

Before LD refinment step, we need to split the genome into smaller regions. Do this by 
running the following command:
```
  cd ${WK_DIR} && python ${PL_DIR}/WEScall.py splitGenome -c user.cfg.yaml
``` 
After running this command, there will generate the folder `${WK_DIR}/phasing` storing the script `${WK_DIR}/phasing/run.sh` and users can submit the genome splitting master job using
```
  cd phasing
  qsub run.sh >> ./logs/submission.log
```

### 6.3. Genotype refinement through phasing

This step performs genotype refinement through phasing by leveraging linkage disequilibrium (LD) information from study samples of external reference panel. Run the following command to generate the job file: 
```
  cd ${WK_DIR} && python ${PL_DIR}/WEScall.py LDRefine -c user.cfg.yaml
```
Before running the phasing procedure, you can modify parameters in the configure file `${WK_DIR}/varCall/cluster.yaml` if necessary. If you want to change the number of segments being phased at one time, please modify the option `N_ARG="--jobs 200"` in line 81 of file `${WK_DIR}/phasing/run.sh`. 

Finally, you can submit the variant phasing master job using
```
  cd phasing & qsub run.sh >> ./logs/submission.log
```
When all above jobs are finished, the genotyping results are stored in `${WK_DIR}/phasing`. For example, users can see genotypes from chromosome 1 in `${WK_DIR}/phasing/1/1.Final.vcf.gz` 

### 6.4. Variant QC

If steps 6.1 through 6.3 have been done successfully, you can perform a series QC procedures described in our paper. Run the following command:
```
  cd ${WK_DIR} && python ${PL_DIR}/WEScall.py QC -c user.cfg.yaml
```
After the QC procedure is finished, the final .vcf files will be located at, e.g., `QC/after_QC/1.after_QC.vcf.gz`. For a list of parameters and their descriptions for the QC procedures, please run `python ${PL_DIR}/WEScall.py QC --help`.

## 7. Frequently used settings and operations

### 7.1. server specific settings
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
### 7.2. Memory settings of variant calling
The joint calling step may take huge memory when the sample size is very large (>1,000). This may kill the program. You can address this issue by either setting the maximum memory usage or splitting genome into smaller regions (default 1Mb).

Increase the maximum memory usage:
```
${PL_DIR}/topMed/scripts/gcconfig.pm Line:27
27  our $batchopts_step2 = "   -l  select=1:ncpus=1:mem=50G -l walltime=24:00:00 -P 13000026 -q production";
```
Split the whole genome into smaller regions (default 1Mb): 
```
${PL_DIR}/topMed/scripts/gcconfig.pm Line:19
19  our $genotypeUnit = 1000000;
```
## 7. Questions
For further questions, pleast raise issues through github (recommended), or contact Jinzhuang Dou <jinzhuangdou198706@gmail.com> or Degang Wu <dwuab@alumni.ust.hk>.
