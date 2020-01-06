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
* Jinzhuang Dou, Degang Wu, Lin Ding, Kai Wang, Minghui Jiang, Xiaoran Chai, Dermot F. Reilly, E Shyong Tai, Jianjun Liu, Xueling Sim, Shanshan Cheng, Chaolong Wang. Using off-target data from whole-exome sequencing to improve genotyping accuracy, association analysis, and phenotype prediction (under review)

## 3. Dependencies
* python (version >= 3.5) with module drmaa (can be installed through `pip`)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (version >= 5.4)
* java (version >= 1.8.0)
* perl (version >= v5.10) with module YAML::XS (can be installed through `cpan`)
* bcftools (version >= 1.9)
* parallel (optional)

## 4. Installation and configuration

You can download our pipeline by the following command:

`git clone https://github.com/dwuab/WEScall.git` 

**Environment variables**

* **PL_DIR**: the directory where the pipeline is located. i.e., the directory of the git-cloned repository. For example, if the cloned repo is at `/opt/software/WEScall`, then `PL_DIR=/opt/software/WEScall`.
* **WK_DIR**: the directory where you run the pipeline and  store the outputs.

### 4.1 Generating 1KG reference panel

**Please run** `${PL_DIR}/scripts/create_g1k_ref.sh` to generate 1000G reference panel files. You should have downloaded 1000G phase 3 data before running this command.

Link to 1000G phase 3 data: [ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)

### 4.2 Downloading resource files

**Please run** `${PL_DIR}/scripts/download_resources.sh` to download resource files ([GotCloud resource bundle](ftp://anonymous@share.sph.umich.edu/gotcloud/ref/hs37d5-db142-v1.tgz) and [ Beagle genetic maps]( http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip)) needed. **This could take a while.** Alternatively, if you already have the resource files downloaded, you can run `${PL_DIR}/scripts/check_resources.sh` to check what resources files you lack and ways to download it. If the script determines a particular resource file is absent, please copy the mentioned resource file to the expected place or make a soft link to it.

### 4.3 Configure the pipeline for your cluster

Review and modify, if necessary, the contents of `${PL_DIR}/cfg/run.template.sh` and `${PL_DIR}/cfg/varCall.cfg.yaml`, according to the cluster engine type, queue name, wall time limits on your cluster.
The default settings are tested on a Torque (An implementation of PBS) cluster. 
We provide two example files, `${PL_DIR}/cfg/run.template.PBSPro.sh` and `${PL_DIR}/cfg/run.template.SGE.sh` to help you set up `${PL_DIR}/cfg/run.template.sh`.
Comments in `${PL_DIR}/cfg/varCall.cfg.yaml` should be helpful too.

## 5. Running the pipeline

### 5.1. Variant calling 
Before running WEScall, you should first prepare a file, `samples.index`, containing a list of samples to call. The file have the same format required by [TopMed](https://github.com/statgen/topmed_freeze3_calling) pipeline.
Each line of `samples.index` is of the following format:

```
  [sampleID] [Absolute Path to BAM/CRAM file] [Contamination rate -- set to zero if unknown].
```
**THERE CANNOT BE EMPTY LINES IN `samples.index` FILE!**
**The path to BAM/CRAM file should be absolute.**
The index file has to be **tab-delimited**. 
BAM/CRAM files listed are assumed to be aligned to `hs37d5.fa`, indexed by `samtools` and **contain no hard clipped reads**, i.e., reads whose CIGAR string contains "H".

You should also prepare a configure file specifying the chromosomes to call, the paths to the resources required by the pipeline and the type of sequencing data (WES or WGS). One example `user.cfg.yaml` in `${PL_DIR}/example/test_WES` is as following: 

```yaml
  chrs: 1,2,10,20,X  
  targetBed:  ${PL_DIR}/resources/SeqCap_EZ_Exome_v3_primary.bed
  1KG3_panel: ${PL_DIR}/resources/1000G_ref_panel/
  geneticMap: ${PL_DIR}/resources/geneticMap_GRCh37/
  seqType:    WES
```
The first line specifies which chromosomes you want to call. Chromosomes that can be called are chromosomes 1 to 22 and X. WEScall can not call variants from chromosome Y. You can specify multiple chromosomes one by one, delimited by comma (for example 20,22). 

The second line specifies the target region bed file. It lists the targeted exonic regions with start and stop chromosome locations in GRCh37/hg19. Note, WEScall can also support the analysis of WGS samples, in which case the target region bed file is not necessary. 

The third line specifies the location of the 1KG3 reference panel used for genotype phasing. See section 4.1.

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

### 5.2. LD-based genotype refinement through phasing

This step performs genotype refinement through phasing by leveraging linkage disequilibrium (LD) information from study samples of external reference panel.
After step 5.1 has done, run the following command to generate the job file: 

```
cd ${WK_DIR} && python ${PL_DIR}/WEScall.py LDRefine -c user.cfg.yaml
```

Finally, you can submit the variant phasing master job using
```
cd LDRefine && qsub run.sh >> ./logs/submission.log
```
When all above jobs are finished, the genotyping results are stored in `${WK_DIR}/LDRefine`. For example, users can see genotypes from chromosome 1 in `${WK_DIR}/LDRefine/1/1.Final.vcf.gz` 

### 5.3. Variant QC

If steps 5.1 and 5.2 have been done successfully, you can perform a series QC procedures described in our paper. Run the following command:
```
  cd ${WK_DIR} && python ${PL_DIR}/WEScall.py QC -c user.cfg.yaml
```
After the QC procedure is finished, the final .vcf files will be located at, e.g., `QC/after_QC/1.after_QC.vcf.gz`. For a list of parameters and their descriptions for the QC procedures, please run `python ${PL_DIR}/WEScall.py QC --help`.

## 6. Frequently used settings and operations

### 6.1. server specific settings
If you want to modify the queue names and pass other parameters to `qsub`. You can modify header of `${PL_DIR}/cfg/run.template.sh`, and `batchopts_step1`, `batchopts_step2`, `batchopts_step3` options of `${PL_DIR}/cfg/varCall.cfg.yaml`.

### 6.2. Memory settings of variant calling
The joint calling step may take huge memory when the sample size is very large (>1,000). The cluster engine may terminate the jobs due to excessive memory usage. You can address this issue by either modifying the amount of requested memory or splitting genome into smaller regions (default 1Mb).

To adjust the maximum memory usage, modify `batchopts_step1`, `batchopts_step2`, `batchopts_step3` options of `${PL_DIR}/cfg/varCall.cfg.yaml`.

To split the whole genome into smaller regions (default 1Mb), modify `genotypeUnit` of `${PL_DIR}/cfg/varCall.cfg.yaml`.

You can also modify the memory and time limits for the jobs in `${PL_DIR}/cfg/cluster.varCall.yaml` and `${PL_DIR}/cfg/cluster.LDRefine.yaml`.

### 6.3 Run the pipeline in local mode
Sometimes it is tricky to set up the pipeline for your own cluster engine, but you want to try this pipeline anyway, you can run the pipeline in local mode.
For `varCall` step, you should set `batchtype` to `local` in `${PL_DIR}/cfg/cluster.varCall.yaml`, and run the step by `./run.sh` instead of `qsub run.sh`.
For `LDRefine` step, simply run `./run.sh`.
`QC` step runs in local mode by default.

## 7. Frequently encountered problems

### 7.1 Pipeline stops prematurely

Symptoms: `qstat` shows the master job as finished, but in the master log you can't find statement `(100%) done`. 
The first thing to do is to find out whether the pipeline has encountered any error in its execution. For example, you can run `grep error ${WK_DIR}/varCall/logs/*` in the log folder to see all mentions of errors. 
See if the error messages come from a particular script or from Snakemake.
See if the error messages clearly point out the underlying sources of errors and if yes try to address the errors.

### 7.2 Network files system synchronization latency

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
but cannot find out other reasons for premature termination, file system latency might be the issue.
Simple resubmission of the job could solve the issue.
Alternatively, you can also increase the time latency by adjusting `time_latency_job` of `${PL_DIR}/cfg/varCall.cfg.yaml`, and modify `DEFAULT_SNAKEMAKE_ARGS` in `${PL_DIR}/cfg/run_template.sh`.

### 7.3 Error message `Can't locate YAML/XS.pm in @INC`

Clearly you don't have perl module YAML::XS installed on your system. Run `cpan install YAML::XL`.

### 7.4 `ModuleNotFoundError: No module named 'drmaa'`

You can install the required python module `drmaa` by `pip install drmaa`.

### 7.5 Mysterious perl error messages

If you get the following message:

```
perl: symbol lookup error: /mnt/software/lib/perl5/5.10.1/auto/Cwd/Cwd.so: undefined symbol: Perl_Istack_sp_ptr
```

it's likely that you have at least two perl installations in your system and one perl installation is trying to load extensions compiled for another perl installation. Set your environment variable `PERL5LIB` appropriately so that versions of the extensions match the perl executable. 

According to our own experience, perl that comes with `miniconda` might trigger such error. If this is the case, try to set `PATH` environment variable so that perl that comes with the Linux distro is used. A workaround is, create a folder `perl_bin` that contains links to `perl` and `cpan`, and then set `PATH` in `~/.bashrc` in the following fashion:

```bash
export PATH=/path/to/perl_bin/:/path/to/conda/:$PATH
```



### 7.6 `/lib64/libc.so.6: version "GLIBC_2.14" not found`

If you encounter error message similar to the one shown above, it's likely that you are running the pipeline on a very old Linux distro, so that some programs such as `samtools`, `bgzip` can't locate the version of the library (or of higher version). We have tested our pipeline on a cluster running CentOS 6.5 with `gcc` 4.4. Distro older than CentOS 6.5 might not be able to run this pipeline.

## 8. Questions
For further questions, please raise issues through github (recommended), or contact Degang Wu <dwuab@alumni.ust.hk> or Jinzhuang Dou <jinzhuangdou198706@gmail.com>.
