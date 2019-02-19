pipelineDir="/mnt/projects/wangcl/ancestry/software/pipelines_web_dev"

# Load the environment variable

pythonPath="/mnt/projects/rpd/apps/miniconda3/bin/python"

# Run BWA


############################################ Run TopMed
# Prepare for the sample configure file
# test.index : aligned sequenced reads in BAM/CRAM format. The third column denotes the contamination rate. 
# test.ped : a pedigree file of nuclear families. As no pedigree information is required in our pipeline, we set the first three columns the sample IDs, the fourth and fifth columns zero. 
# user.cfg.yaml : specify parameters 

rm -rf varCall

# generate the command 
${pythonPath} ${pipelineDir}/varCall/runTopMed.py -c  user.cfg.yaml   -i test.index -p test.ped  -t WGS -n
cd varCall & qsub run.sh >> ./logs/submission.log

######################################### Split the genome into small regions 
rm -rf phasing

${pythonPath} ${pipelineDir}/phasing/splitGenome.py  -t WGS -c user.cfg.yaml -n
cd phasing & qsub run.sh >> ./logs/submission.log

######################################### Run Beagle 
${pythonPath} ${pipelineDir}/phasing/phasing.py  -t WGS -c user.cfg.yaml -n
cd phasing & qsub run.sh >> ./logs/submission.log

