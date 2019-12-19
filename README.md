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
* Jinzhuang Dou, Degang Wu, Lin Ding, Kai Wang, Minghui Jiang, Xiaoran Chai, Dermot F. Reilly, E Shyong Tai, Jianjun Liu, Xueling Sim, Shanshan Cheng, Chaolong Wang. Using off-target data from whole-exome sequencing to improve genotyping accuracy, association analysis, and phenotype prediction (in preparation)

## 3. Dependencies
* python (version >= 3.5)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (version >= 5.4)
* java (version >= 1.8.0)
* perl (version >= v5.10)
* perl module YAML::XS
* bcftools (version >= 1.9)
* parallel (optional)

## 4. Download from github

You can download our pipeline by the following command:

`git clone https://github.com/dwuab/WEScall.git` 

## 5. Getting started 

**Environment variables**
* **PL_DIR**: the directory where the pipeline is located. 
* **WK_DIR**: the directory where you store the outputs.

### 5.1 Setting up for your cluster

Review and modify, if necessary, the contents of `${PL_DIR}/cfg/run.template.sh` and `${PL_DIR}/cfg/varCall.cfg.yaml`, according to the cluster engine type, queue name, wall time limits on your cluster.
The default settings are tested on a Torque (An implementation of PBS) cluster. 
We provide two example files, `${PL_DIR}/cfg/run.template.PBSPro.sh` and `${PL_DIR}/cfg/run.template.SGE.sh` to help you set up `${PL_DIR}/cfg/run.template.sh`.
Comments in `${PL_DIR}/cfg/varCall.cfg.yaml` should be helpful too.

### 5.2 Generating 1KG reference panel

**Please run** `${PL_DIR}/scripts/create_g1k_ref.sh` to generate 1000G reference panel files. You should have downloaded 1000G phase 3 data before running this command.

### 5.3 Setting up links to resource files

**Please run** `${PL_DIR}/scripts/check_resources.sh` to check what resources files you lack and ways to download it. If the script determines the absence of a particular resource file, please copy the mentioned resource file to the expected place or make a soft link to it.

## 6. Running the pipeline

### 6.1. Variant calling 
Before running WEScall, you should first prepare a file, `samples.index`, containing a list of samples to call. The file have the same format required by [TopMed](https://github.com/statgen/topmed_freeze3_calling) pipeline.
Each line of `samples.index` is of the following format:
```
  [sampleID] [Absolute Path to BAM/CRAM file] [Contamination rate -- set to zero if unknown].
``` 
**THERE CANNOT BE EMPTY LINES IN `samples.index` FILE!**
**The path to BAM/CRAM file should be absolute.**
The index file has to be **tab-delimited**. 
BAM/CRAM files listed are assumed to be indexed and **contain no hard clipped reads**, i.e., reads whose CIGAR string contains "H".

You should also prepare a configure file specifying the chromosomes to call, the paths to the resources required by the pipeline and the type of sequencing data (WES or WGS). One example `user.cfg.yaml` in `${PL_DIR}/example/test_WES` is as following: 

```
  chrs: 1,2,10,20,X  
  targetBed:  ${PL_DIR}/WEScall/resources/SeqCap_EZ_Exome_v3_primary.bed
  1KG3_panel: ${PL_DIR}/WEScall/resources/1000G_ref_panel
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
cd varCall && qsub run.sh >> ./logs/submission.log  
``` 
If any job is killed prematurely, you can resume the master job by using the command again. You can check `${WK_DIR}/varCall/logs/WEScall_varCall.master.log` for progress or diagnose premature terminations of jobs. 
If the running is successful, the vcf files after SVM filtering are placed in, e.g.,  `${WK_DIR}/varCall/1/1.Filter.vcf.gz`

Once the log reports all jobs are done (message such as "4 of 4 steps (100%) done"), you can proceed to the next step.

### 6.2. LD-based genotype refinement through phasing

This step performs genotype refinement through phasing by leveraging linkage disequilibrium (LD) information from study samples of external reference panel.
After step 6.1 has done, run the following command to generate the job file: 
```
cd ${WK_DIR} && python ${PL_DIR}/WEScall.py LDRefine -c user.cfg.yaml
```

Finally, you can submit the variant phasing master job using
```
cd phasing && qsub run.sh >> ./logs/submission.log
```
When all above jobs are finished, the genotyping results are stored in `${WK_DIR}/LDRefine`. For example, users can see genotypes from chromosome 1 in `${WK_DIR}/LDRefine/1/1.Final.vcf.gz` 

### 6.3. Variant QC

If steps 6.1 and 6.2 have been done successfully, you can perform a series QC procedures described in our paper. Run the following command:
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

## 8. Frequently encountered problems

### 8.1 Pipeline stops prematurely

Symptoms: `qstat` shows the master job as finished, but in the master log you can't find statement `(100%) done`. 
The first thing to do is to find out whether the pipeline has encountered any error in its execution. For example, you can run `grep error ${WK_DIR}/varCall/logs/*` in the log folder to see all mentions of errors. 
See if the error messages come from a particular script or from Snakemake.
See if the error messages clearly point out the underlying sources of errors and if yes try to address the errors.

### 8.2 Network files system synchronization latency

Occasionally, a job has finished, but it takes a long time for the outputs it generated to be synchronized to other computing nodes. 
In this case, the master job will be informed by the cluster scheduler that the job has been finished but unable to detect the expected output files.
For example, if you see sentences like the following in the log file:
```
Waiting at most 600 seconds for missing files.
MissingOutputException in line 191 of /opt/software/WEScall/pipelines/LDRefine/Snakefile.beagle.WES:
Missing files after 600 seconds:
20_split/20
This might be due to filesystem latency.
```
but cannot find out the 

## 9. Questions
For further questions, please raise issues through github (recommended), or contact Degang Wu <dwuab@alumni.ust.hk> or Jinzhuang Dou <jinzhuangdou198706@gmail.com>.
