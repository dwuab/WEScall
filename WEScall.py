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

from pipelines import get_cluster_cfgfile
from pipelines import PipelineHandler

# global logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
	'[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)


def get_seq_type_from_user_cfg(fn_user_cfg):
	with open(fn_user_cfg) as fh:
		user_cfg = dict(yaml.safe_load(fh))
	seq_type=user_cfg["seqType"]
	if (seq_type!="WES") and (seq_type!="WGS"):
		logger.error("seqType can only be WES or WGS! exiting!")
		exit(1)

	return seq_type

def varCall(args):
	logger.debug("Varcall: "+str(args))

	PIPELINE_BASEDIR = os.path.dirname(sys.argv[0])
	CFG_DIR = os.path.join(PIPELINE_BASEDIR, "cfg")

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
	logger.debug("LDRefine: "+str(args))

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
	logger.debug("QC: "+str(args))

	assert (args.DR2>=0. and args.DR2<=1.0), "DR2 has to be a floating point number between 0 and 1!"

	assert args.MAC>=0, "MAC has to be an integer larger than or eqaul to 0."

	assert (args.HWE>=0. and args.HWE<=1.), "HWE p-value has to be between 0 and 1."

	assert args.DP_target>0., "DP_target has to be larger than 0."

	assert args.DP_off_target>0., "DP_off_target has to be larger than 0."

	assert args.indel_excl>0, "indel_excl has to be an integer larger than or equal to 0."

	assert (args.max_GP>=0. and args.max_GP<=1.0), "max_GP has to be between 0 and 1."

	assert (args.missing_rate>=0. and args.missing_rate<=1.0), "missing_rate has to be between 0 and 1."

	if get_seq_type_from_user_cfg(args.userCfg)=="WGS":
		logger.warn("This QC procedure was designed for WES data.")

	if not os.path.exists("QC"):
		os.mkdir("QC")

	with open(args.userCfg) as fh_cfg:
		user_cfg_dict=yaml.safe_load(fh_cfg)

	user_cfg_dict["QC"]={"DR2":args.DR2,
	"MAC":args.MAC,
	"HWE":args.HWE,
	"DP_target":args.DP_target,
	"DP_off_target":args.DP_off_target,
	"g1k_indel_interval_bed":args.g1k_intdel_interval_bed,
	"indel_excl":args.indel_excl,
	"max_GP":args.max_GP,
	"missing_rate":args.missing_rate
	}

	with open("./QC/QC_params.yaml","w") as fh_cfg:
		fh_cfg.write(yaml.dump(user_cfg_dict))

	PIPELINE_BASEDIR = os.path.dirname(sys.argv[0])

	cmd="PL_DIR="+PIPELINE_BASEDIR+" && cd QC && snakemake -s ${PL_DIR}/pipelines/QC/Snakefile.QC.WES "+\
		"--configfile QC_params.yaml --cores 20 --printshellcmds all"
	logger.debug(cmd)
	os.system(cmd)


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

	# one subparser corresponds to one subcommand
	parser_varCall = subparsers.add_parser('varCall', parents=[common_parser],
		help='Variant discovery, genotype calling and site filtering by SVM')
	parser_varCall.add_argument('-s', '--sample-list', required=True,
								help="The list of study samples, paths to the corresponding BAM/CRAM files and contamination rates")
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
	parser_QC.add_argument('--DP_target',type=float,default=300,metavar="depth",
		help="Depth threshold for target region. Range: integer>=0. Variants in the target region with DP larger than this value will be filtered.")
	parser_QC.add_argument('--DP_off_target',type=float,default=50,metavar="depth",
		help="Depth threshold for target region. Range: integer>=0. Variants in the off-target region with DP larger than this value will be filtered.")
	parser_QC.add_argument('--g1k_intdel_interval_bed',type=str,metavar="path",
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
	assert os.path.exists(args.userCfg), ("User config file {} cannot be found".format(args.userCfg))

	# execute subcommand-specific function
	args.func(args)

if __name__ == "__main__":
	main()
