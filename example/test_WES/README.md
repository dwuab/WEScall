This file guides you to run the WES example

step 1: download the public WES data from Genome In A Bottle (GIAB) project. If you already have this data, you can make soft links to them in this folder.

run: 
```
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bam
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup.bam
```
samtools index 151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam

step 2: index the downloaded bam files
```
samtools index 151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
samtools index 151002_7001448_0359_AC7F6GANXX_Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bam
samtools index 151002_7001448_0359_AC7F6GANXX_Sample_HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup.bam
```

step 3: set up the pipeline, which is already described in section 4 of the main README.md

step 4: make sure the dependencies of the pipeline have been satisfied. ***In particular, make sure python3 can be found in $PATH.***

step 5: review and modify parameters, **in particular, paths to various resource files**, in `user.cfg.yaml` if necessary. Do the same to `sample.index`, and in particular, use absolute paths for all bam files. Set the environment variable `PL_DIR=/path/to/WEScall`.

step 6: run `python ${PL_DIR}/WEScall.py varCall -c user.cfg.yaml -s samples.index`. If no error is encountered. Run `cd varCall && qsub run.sh >> ./logs/submission.log`. Check the status of jobs through `qstat` or equivalent commands. If you find the phrase `100% done` at the end of `varCall/logs/WEScall_varCall.master.log`. This step is done.

step 7: run `python ${PL_DIR}/WEScall.py LDRefine -c user.cfg.yaml`. The rest is similar to the previous step. Check `LDRefine/logs/WEScall_LDRefine.master.log` for the status of the job.

step 8: run `python ${PL_DIR}/WEScall.py QC -c user.cfg.yaml` until everything is done. Review the final result at `QC/after_QC/20.after_QC.vcf.gz`.

step 9: You can use `bcftools stats 20.after_QC.ref.vcf.gz QC/after_QC/20.after_QC.vcf.gz -s HG002,HG003,HG004` to inspect the non-reference discordance rate.
