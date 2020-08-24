#!/usr/bin/env python3
"""
The main interface of WEScall
"""

import argparse
import sys
import os
import yaml
import logging
import shutil
import glob
import re

LIB_PATH = os.path.abspath(
	os.path.join(os.path.dirname(os.path.realpath(__file__)), "pipelines/lib"))

if LIB_PATH not in sys.path:
	sys.path.insert(0, LIB_PATH)

PIPELINE_BASEDIR = os.path.dirname(os.path.realpath(sys.argv[0]))
CFG_DIR = os.path.join(PIPELINE_BASEDIR, "cfg")

import pipelines
from pipelines import get_cluster_cfgfile
from pipelines import PipelineHandler

# global logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
	'[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)


def print_parameters_given(args):
	logger.info("Parameters in effect:")
	for arg in vars(args):
		if arg=="func": continue
		logger.info("--{} = [{}]".format(arg, vars(args)[arg]))


def get_seq_type_from_user_cfg(fn_user_cfg):
	with open(fn_user_cfg) as fh:
		user_cfg = dict(yaml.safe_load(fh))
	seq_type=user_cfg["seqType"]
	if (seq_type!="WES") and (seq_type!="WGS"):
		logger.error("seqType can only be WES or WGS! exiting!")
		exit(1)

	return seq_type


def validate_sample_list_file(args):
	if args.check_hard_clipped:
		out=os.popen("command -v bioawk").read().strip()
		assert out!="", "Program bioawk cannot be found!"

	assert os.path.isfile(args.sample_list), "Sample index file {} cannot be found!".format(args.sample_list)

	try:
		with open(args.sample_list) as f_in:
			for line in f_in:
				record = line.strip().split("\t")
				logger.debug("Checking sample {}".format(record[0]))
				assert len(record)==3, "Every line has to have exactly 3 tab-delimited columns! Line with sample name {} does not satisify this requiremnt!".format(record[0])
				assert os.path.isfile(record[1]), "Bam file {} cannot be found!".format(record[1])
				assert os.path.isfile(record[1]+".bai"), "Bam file {} has not been indexed!".format(record[1])
				assert os.path.isabs(record[1]), "Please use absolute path for bam file {}!".format(record[1])

				if args.check_hard_clipped:
					logger.debug("Checking existence of hard-clipped reads.")
					cmd = "samtools view {} | bioawk -c sam 'BEGIN {{count=0}} ($cigar ~ /H/)&&(!and($flag,256)) {{count++}} END {{print count}}'".format(record[1])
					logger.debug("Command: "+cmd)
					out=os.popen(cmd).read().strip()
					logger.debug("Results: "+out)
					assert out=="0", "Bam file {} contains hard-clipped reads without proper flag (0x100) set! Please use -M or -Y options of BWA MEM!".format(record[1])

				try:
					float(record[2])
					assert 0.0 <= float(record[2]) and float(record[2]) <= 1.0, "Contamination rate of sample {0} has to be a float number between 0 and 1 instead of {1}!".format(record[0], record[2])
				except:
					logger.error("Contamination rate of sample {0} has to be a float number between 0 and 1 instead of {1}!".format(record[0], record[2]))
					exit(1)

	except Exception:
		logger.error("There is something wrong with the sample index file. Check the logs for more information.")
		print(sys.exc_info())
		raise sys.exc_info()[0]


def validate_user_cfg(args):
	with open(args.userCfg) as fh:
		user_cfg = dict(yaml.safe_load(fh))

	if user_cfg["seqType"]=="WES":
		assert os.path.isfile(user_cfg["targetBed"]), "Target region bed file {} cannot be found!".format(user_cfg["targetBed"])

	if user_cfg["seqType"]=="WGS":
		chrs = [item.strip() for item in str(user_cfg["chrs"]).split(",")]
		assert ("1" in chrs) and ("X" in chrs), "Chromosomes 1 and X have to be included in WGS mode!"


	assert os.path.isdir(user_cfg["1KG3_panel"]), "1KG3 reference panel cannot be found at {}!".format(user_cfg["1KG3_panel"])
	all_chr=list(range(1,23))
	all_chr.extend("X")
	for chr_no in all_chr:
		site_fn=os.path.join(user_cfg["1KG3_panel"],"ALL.chr{}.phase3.20130502.SNP.biallelic.MAF0.01.sites.vcf.gz".format(chr_no))
		site_index_fn=os.path.join(user_cfg["1KG3_panel"],"ALL.chr{}.phase3.20130502.SNP.biallelic.MAF0.01.sites.vcf.gz.csi".format(chr_no))
		GT_fn=os.path.join(user_cfg["1KG3_panel"],"ALL.chr{}.phase3.20130502.SNP.indel.biallelic.mac5.vcf.gz".format(chr_no))
		
		assert os.path.isfile(site_fn), "1KG3 ref panel site file {} cannot be found!".format(site_fn)
		assert os.path.isfile(site_index_fn), "Index file of 1KG3 ref panel site file {} cannot be found!".format(site_fn)
		assert os.path.isfile(GT_fn), "Filtered genotype file of 1KG3 ref panel {} cannot be found!".format(GT_fn)

	assert os.path.isdir(user_cfg["geneticMap"]), "Genetic map folder cannot be found at {}.".format(user_cfg["geneticMap"])

	for chr_no in all_chr:
		map_fn=os.path.join(user_cfg["geneticMap"], "plink.chr{}.GRCh37.map".format(chr_no))
		assert os.path.isfile(map_fn), "Genetic map file {} cannot be found!".format(map_fn)


def check_resource_files_for_varCall():
	logger.debug(PIPELINE_BASEDIR)
	resource_folder=os.path.join(PIPELINE_BASEDIR, "pipelines/varCall/gotcloud.ref")

	assert os.path.isdir(resource_folder), "Resource folder {} cannot be found!".format(resource_folder)

	required_resource_files=("hs37d5.fa.fai", 
		"hs37d5.fa", 
		"1000G_omni2.5.b37.sites.PASS.vcf.gz",
		"1000G_omni2.5.b37.sites.PASS.vcf.gz.tbi",
		"hapmap_3.3.b37.sites.vcf.gz",
		"hapmap_3.3.b37.sites.vcf.gz.tbi",
		"dbsnp_142.b37.vcf.gz",
		"dbsnp_142.b37.vcf.gz.tbi")

	for fn in required_resource_files:
		fn_path=os.path.join(PIPELINE_BASEDIR, "pipelines/varCall/gotcloud.ref", fn)
		assert os.path.isfile(fn_path), "Required resource file {} cannot be found!".format(fn_path)


def check_dependencies():
	programs_to_check = ("snakemake", "perl", "bcftools", "java")

	for prog in programs_to_check:
		out = os.popen("command -v {}".format(prog)).read()
		assert out != "", "Program {} cannot be found!".format(prog)

	perl_mods_to_check = ("YAML::XS",)

	for mod in perl_mods_to_check:
		out_pipe = os.popen('perl -e "use {}"'.format(mod))

		# if the specified perl module has been installed, then there will be no output,
		# otherwise out.close() will returns an error number

		assert out_pipe.close() is None, "Perl module {} has not been installed!".format(mod)

#	python_pkgs_to_check = ("drmaa",)

#	for pkg in python_pkgs_to_check:
#		out_pipe = os.popen('python -c "import {}"'.format(pkg))

#		assert out_pipe.close() is None, "Python module {} has not been installed!".format(pkg)


def varCall(args):
	logger.info("Preparing varCall pipeline...")
	print_parameters_given(args)

	logger.info("Validating sample index ...")
	validate_sample_list_file(args)

	logger.info("Validating user config file ...")
	validate_user_cfg(args)

	logger.info("Checking existence of essenstial resource files...")
	check_resource_files_for_varCall()

	logger.info("Checking dependencies...")
	check_dependencies()

	pipeline_handler = PipelineHandler(
		"WEScall_varCall",
		PIPELINE_BASEDIR,
		Snakefile="pipelines/varCall/Snakefile."+get_seq_type_from_user_cfg(args.userCfg),
		outdir="./varCall",
		user_data="",
		user_cfgfile=args.userCfg,
		cluster_cfgfile=CFG_DIR+"/cluster.varCall.yaml"
		)


	os.system("mkdir -p ./varCall/data")
	shutil.copy2(args.sample_list,"./varCall/data/samples.index")

	# automatically generate the pedigree file for the user
	# Since WEScall does not utilize pedigree information, the pedigree file
	# is just a formality so as to let the pipeline run
	with open(args.sample_list) as f_in, open("./varCall/data/samples.ped","w") as f_out:
		for line in f_in:
			record = line.strip().split("\t")
			f_out.write("{smp}\t{smp}\t{smp}\t0\t0\n".format(smp=record[0]))

	pipeline_handler.setup_env()
	pipeline_handler.submit(no_run=True)


def LDRefine(args):
	logger.info("Preparing LD-based genotype refinement pipeline...")
	print_parameters_given(args)

	assert os.path.exists("varCall"), "Cannot detect the directory of varaiant detection.\nWEScall varCall has to be run before LD-based genotype refinement."

	assert args.num_record_per_file > 0, "Number of records per file has to be larger than 1!"

	assert args.num_overlap_record >= 0, "Number of overlapping records has to be larger than 1!"

	assert args.num_record_per_file>args.num_overlap_record, "Number of records per file has to be larger than the number of overlapping records."

	if not os.path.exists("LDRefine"):
		os.mkdir("LDRefine")

	LDRefine_cfg=dict()
	LDRefine_cfg["num_record_per_file"]=args.num_record_per_file
	LDRefine_cfg["num_overlap_record"]=args.num_overlap_record

	PIPELINE_BASEDIR = os.path.join(os.path.dirname(sys.argv[0]))
	CFG_DIR = os.path.join(PIPELINE_BASEDIR, "cfg")

	path_cluster_cfg=os.path.join(PIPELINE_BASEDIR,"cfg","cluster.LDRefine.yaml")

	# has to merge cluster
	with open(path_cluster_cfg, 'r') as fh:
		cluster_cfg = yaml.safe_load(fh)

	# turn arguments into user_data that gets merged into pipeline config
	#
	# generic data first
	user_data = dict()
	user_data['cluster'] = cluster_cfg
	user_data['LDRefine'] = LDRefine_cfg

	pipeline_handler = PipelineHandler(
		"WEScall_LDRefine", PIPELINE_BASEDIR, 
		"LDRefine",user_data,
		Snakefile="pipelines/LDRefine/Snakefile.beagle."+get_seq_type_from_user_cfg(args.userCfg),
		cluster_cfgfile=path_cluster_cfg,
		user_cfgfile=args.userCfg)

	pipeline_handler.setup_env()
	pipeline_handler.submit(no_run=True)

def QC(args):
	logger.info("Preparing final QC pipeline...")
	print_parameters_given(args)

	assert (args.DR2>=0. and args.DR2<=1.0), "DR2 has to be a floating point number between 0 and 1!"

	assert args.MAC>=0, "MAC has to be an integer larger than or eqaul to 0."

	assert (args.HWE>=0. and args.HWE<=1.), "HWE p-value has to be between 0 and 1."

	assert (args.DP_target_upper>0. and args.DP_target_lower>=0. and args.DP_target_upper>args.DP_target_lower), "DP_target_lower and DP_target_upper have to be larger than 0, and DP_target_upper has to be larger than DP_target_lower."

	assert args.DP_off_target>0., "DP_off_target has to be larger than 0."

	assert args.indel_excl>0, "indel_excl has to be an integer larger than or equal to 0."

	assert (args.max_GP>=0. and args.max_GP<=1.0), "max_GP has to be between 0 and 1."

	assert (args.missing_rate>=0. and args.missing_rate<=1.0), "missing_rate has to be between 0 and 1."

	if ("female_sample_list" in args):
		assert (os.path.isfile(args.female_sample_list)), "female sample list {} cannot be found!".format(args.female_sample_list)
		if (args.skip_HWE_X):
			logger.error("A female sample list is supplied while user choose to skip HWE filtering on X chromosome. "
				"Either set --skip_HWE_X and not --female_sample_list, or not set --skip_HWE_X and set --female_sample_list."
				"... Exiting")
			exit()

	if get_seq_type_from_user_cfg(args.userCfg)=="WGS":
		logger.warn("This QC procedure was designed for WES data.")

	if not os.path.exists("QC"):
		os.mkdir("QC")

	with open(args.userCfg) as fh_cfg:
		user_cfg_dict=yaml.safe_load(fh_cfg)

	user_cfg_dict["QC"]={
	"DR2":args.DR2,
	"MAC":args.MAC,
	"HWE":args.HWE,
	"DP_target_upper":args.DP_target_upper,
	"DP_target_lower":args.DP_target_lower,
	"DP_off_target":args.DP_off_target,
	"g1k_indel_interval_bed":args.g1k_indel_interval_bed,
	"indel_excl":args.indel_excl,
	"max_GP":args.max_GP,
	"missing_rate":args.missing_rate,
	"skip_HWE_X":args.skip_HWE_X
	}

	if ("female_sample_list" in args):
		user_cfg_dict["QC"]["female_sample_list"] = args.female_sample_list

	with open("./QC/QC_params.yaml","w") as fh_cfg:
		fh_cfg.write(yaml.dump(user_cfg_dict))

	PIPELINE_BASEDIR = os.path.dirname(os.path.realpath(__file__))

	cmd="set -euo pipefail && PL_DIR="+PIPELINE_BASEDIR+" && cd QC && snakemake -s ${PL_DIR}/pipelines/QC/Snakefile.QC.WES "+\
		"--configfile QC_params.yaml --cores 20 --printshellcmds all 2>&1 | tee QC.log"
	logger.debug(cmd)
	logger.info("Start running QC procedures.")
	exit_code=os.system(cmd)
	if (exit_code!=0):
		logger.error("Errors encountered. Exit code: {}. See logs above for information.".format(exit_code))
		exit(exit_code)
	

def main():
	parser = argparse.ArgumentParser(
		description="""A whole exome sequence genotype calling pipeline that utilizes off-target reads)
		""",
		epilog=
		"""Typical workflow: varCall => LDRefine => QC\nCitation: Dou J, Wu D, ... , Wang C. Joint analysis of target and off-target data improves whole-exome sequencing studies (in preparation)
		""",
		formatter_class=argparse.RawTextHelpFormatter)

	subparsers = parser.add_subparsers(title='Available subcommands', dest="subcommand")
	
	# every subcommand needs user config file
	common_parser = argparse.ArgumentParser(add_help=False)
	common_parser.add_argument('-c', '--userCfg', required=True,
						metavar="userCfgFile",
						help="User specific configure file")
	common_parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")

	# one subparser corresponds to one subcommand
	parser_varCall = subparsers.add_parser('varCall', parents=[common_parser],
		help='Variant discovery, genotype calling and site filtering by SVM',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_varCall.add_argument('-s', '--sample-list', required=True,
								help="The list of study samples, paths to the corresponding BAM/CRAM files and contamination rates")
	parser_varCall.add_argument('-H', '--check-hard-clipped',
		help="See existence of hard-clipped reads with improper flag (0x100 not set) in bam files. Requires program bioawk.",
		action="store_true")
	parser_varCall.set_defaults(func=varCall)

	parser_LDRefine = subparsers.add_parser('LDRefine', parents=[common_parser],
		help='Perform LD-based genotype refinement',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_LDRefine.add_argument("--num-record-per-file",type=int,default=10000,
		metavar="number-of-vcf-records-per-file",
		help="Number of VCF records per split file. Range: integer>0.")
	parser_LDRefine.add_argument("--num-overlap-record",type=int,default=1000,
		metavar="number-of-overlapping-vcf-records",
		help="Number of overlapping VCF records between two neighboring split vcf files. Range: integer>=0.")
	parser_LDRefine.set_defaults(func=LDRefine)
	
	parser_QC = subparsers.add_parser('QC', parents=[common_parser],
		help='final quality control',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_QC.add_argument('--DR2',type=float,default=0.8,metavar="float",
		help="DR2 threshold. Range: [0,1]. Variants with DR2 less than this value will be filtered.")
	parser_QC.add_argument('--MAC',type=int,default=0,metavar="integer",
		help="MAC threshold. Range: integer>=0. Variants with MAC less than or equal to this value will be filtered.")
	parser_QC.add_argument('--HWE',type=float,default=1e-5,metavar="HWE p-value",
		help="HWE p value threshold. Range: [0,1]. Variants with HWE p-value less than this value will be filtered.")
	parser_QC.add_argument('--skip_HWE_X',action="store_true",default=False,
		help="Whether to skip HWE filtering on X chromosome.")
	parser_QC.add_argument('--female_sample_list',type=str,metavar="file name",
		help="List of female samples for HWE filtering of X chromosome.")
	
	parser_QC.add_argument('--DP_target_upper',type=float,default=300,metavar="depth",
		help="Depth threshold for target region. Range: integer>=0. Variants in the target region with DP larger than this value will be filtered.")
	parser_QC.add_argument('--DP_target_lower',type=float,default=0,metavar="depth",
		help="Depth threshold for target region. Range: integer>=0. Variants in the target region with DP smaller than this value will be filtered.")
	parser_QC.add_argument('--DP_off_target',type=float,default=50,metavar="depth",
		help="Depth threshold for target region. Range: integer>=0. Variants in the off-target region with DP larger than this value will be filtered.")
	
	parser_QC.add_argument('--g1k_indel_interval_bed',type=str,metavar="path",
		default=os.path.join(os.path.dirname(os.path.realpath(__file__)), "resources","g1k.indel.interval.bed.gz"),
		help="Path to 1000G interval list. No need to set this.")
	parser_QC.add_argument('--indel_excl',type=int,default=5,metavar="number-of-bases",
		help="Indel exclusion radius (in number of bases). Range: integer>=0. SNPs within the exclusion radius of indels in 1000G will be filtered.")
	parser_QC.add_argument('--max_GP',type=float,default=0.9,
		help="Maximum GP threshold. Range: [0,1]. If the maximum genotype probability of a sample at a site is less than this value, then this genotype will be considered as missing.")
	parser_QC.add_argument('--missing_rate',type=float,default=0.05,
		help="Missing rate threshold. Range: [0,1]. If the missing rate of a site is less than this value, the site will be filtered.")
	parser_QC.set_defaults(func=QC)

	args = parser.parse_args()
	if args.subcommand is None:
		# if no command is specified, print help and exit
		print("Please specify one subcommand! Exiting!")
		print("-"*80)
		parser.print_help()
		exit(1)

	# check if user configure exists
	assert os.path.isfile(args.userCfg), ("User config file {} cannot be found".format(args.userCfg))

	if args.verbose:
		for name in logging.root.manager.loggerDict:
			logging.getLogger(name).setLevel(logging.DEBUG)

	# execute subcommand-specific function
	args.func(args)

	logger.info("Success! See instructions above.")

if __name__ == "__main__":
	main()
