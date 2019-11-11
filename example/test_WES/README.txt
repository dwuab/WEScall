############################################
# set the paths appropriately according to the settings on your server
pipelineDir=`cwd`/..
pythonPath="/opt/software/miniconda3/bin/python"

############################################ 
# Prepare for the sample configure file
# samples.index : aligned sequenced reads in BAM/CRAM format. The third column denotes the contamination rate. 
# user.cfg.yaml : specify parameters 

# generate the command 
${pythonPath} ${pipelineDir}/varCall/runTopMed.py -c  user.cfg.yaml -i samples.index -t WES -n
cd varCall
qsub run.sh >> ./logs/submission.log

######################################### 
# Split the genome into small regions 

rm -rf phasing

${pythonPath} ${pipelineDir}/phasing/splitGenome.py  -t WES -c user.cfg.yaml -n
cd phasing
qsub run.sh >> ./logs/submission.log

######################################### 
# Run phasing 
${pythonPath} ${pipelineDir}/phasing/phasing.py  -t WES -c user.cfg.yaml -n
cd phasing
qsub run.sh >> ./logs/submission.log

