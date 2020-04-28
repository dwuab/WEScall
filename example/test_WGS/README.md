This file guides you to run the WGS example

step 1: download the public WGS data from Genome In A Bottle (GIAB) project. (These files are large! You may want to download part of the files) If you already have this data, you can make soft links to them in this folder.

run: 

```bash
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.hs37d5.2x250.bam
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/NIST_Illumina_2x250bps/novoalign_bams/HG003.hs37d5.2x250.bam
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/NIST_Illumina_2x250bps/novoalign_bams/HG004.hs37d5.2x250.bam
```

Step 2: index the downloaded bam files

```bash
samtools index HG002.hs37d5.2x250.bam
samtools index HG003.hs37d5.2x250.bam
samtools index HG004.hs37d5.2x250.bam
```

step 3: set up the pipeline, which is already described in section 4 of the main README.md

step 4: review and modify parameters, **in particular, paths to various resource files**, in `user.cfg.yaml` if necessary. **Note that chromosomes 1 and X have to be included.** Do the same to `sample.index`, and in particular, use absolute paths for all bam files. Set the environment variable `PL_DIR=/path/to/WEScall`.

step 5: run `python ${PL_DIR}/WEScall.py varCall -c user.cfg.yaml -s samples.index`. If no error is encountered. Run `cd varCall && qsub run.sh >> ./logs/submission.log`. Check the status of jobs through `qstat` or equivalent commands. If you find the phrase `100% done` at the end of `varCall/logs/WEScall_varCall.master.log`. This step is done.

step 6: run `python ${PL_DIR}/WEScall.py LDRefine -c user.cfg.yaml`. The rest is similar to the previous step. Check `LDRefine/logs/WEScall_LDRefine.master.log` for the status of the job.

step 7: run `python ${PL_DIR}/WEScall.py QC -c user.cfg.yaml` until everything is done. Review the final result at `QC/after_QC/20.after_QC.vcf.gz`.

step 8: You can use `bcftools stats 20.after_QC.ref.vcf.gz QC/after_QC/20.after_QC.vcf.gz -s HG002,HG003,HG004` to inspect the non-reference discordance rate.