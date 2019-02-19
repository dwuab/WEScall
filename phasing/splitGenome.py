#!/usr/bin/env python3
"""{PIPELINE_NAME} pipeline (version: {PIPELINE_VERSION}): creates
pipeline-specific config files to given output directory and runs the
pipeline (unless otherwise requested).
"""
# generic usage {PIPELINE_NAME} and {PIPELINE_VERSION} replaced while
# printing usage

#--- standard library imports
#
import sys
import os
import argparse
import logging

#--- third-party imports
#
import yaml

#--- project specific imports
#
# add lib dir for this pipeline installation to PYTHONPATH
LIB_PATH = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)),"..", "lib"))

if LIB_PATH not in sys.path:
    sys.path.insert(0, LIB_PATH)


from pipelines import get_pipeline_version
from pipelines import PipelineHandler
from pipelines import logger as aux_logger
from pipelines import get_cluster_cfgfile
from pipelines import get_default_queue

__author__ = "Jinzhuang Dou"
__email__ = "douj@gis.a-star.edu.sg"
__copyright__ = "2016 Genome Institute of Singapore"
__license__ = "The MIT License (MIT)"


# only dump() and following do not automatically create aliases
yaml.Dumper.ignore_aliases = lambda *args: True


PIPELINE_BASEDIR = os.path.dirname(sys.argv[0])
CFG_DIR = os.path.join(PIPELINE_BASEDIR, "cfg")

# same as folder name. also used for cluster job names
PIPELINE_NAME = "phasing"

MARK_DUPS = True

# global logger
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)

def main():
    """main function
    """

    parser = argparse.ArgumentParser(description=__doc__.format(
        PIPELINE_NAME=PIPELINE_NAME, PIPELINE_VERSION=get_pipeline_version()))

    # generic args
    parser.add_argument('-c', '--userCfg', required=True, default=0,
                        help="User specific configure file")
    parser.add_argument('--name',
                        help="Give this analysis run a name (used in email and report)")
    parser.add_argument('--no-mail', action='store_true',
                        help="Don't send mail on completion")
    default = get_default_queue('slave')
    parser.add_argument('-w', '--slave-q', default=default,
                        help="Queue to use for slave jobs (default: {})".format(default))
    default = get_default_queue('master')
    parser.add_argument('-m', '--master-q', default=default,
                        help="Queue to use for master job (default: {})".format(default))
    parser.add_argument('-n', '--no-run', action='store_true')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="Increase verbosity")
    parser.add_argument('-q', '--quiet', action='count', default=0,
                        help="Decrease verbosity")
    cfg_group = parser.add_argument_group('Configuration files (advanced)')
    #cfg_group.add_argument('--sample-cfg',
    #                       help="Config-file (YAML) listing samples and readunits."
    #                       " Collides with -1, -2 and -s")
    for name, descr in [("references", "reference sequences"),
                        ("params", "parameters"),
                        ("modules", "modules")]:
        default = os.path.abspath(os.path.join(CFG_DIR, "{}.yaml".format(name)))
        cfg_group.add_argument('--{}-cfg'.format(name),
                               default=default,
                               help="Config-file (yaml) for {}. (default: {})".format(descr, default))

    # pipeline specific args
    parser.add_argument('-t', "--seqtype", required=True,
                        choices=['WGS', 'WES'],
                        help="Sequencing type")
    parser.add_argument('-l', "--intervals",
                        help="Intervals file (e.g. bed file) listing regions of interest."
                        " Required for WES and targeted sequencing.")

    args = parser.parse_args()

    default = os.path.abspath(args.userCfg)

    cfg_group.add_argument('--{}-cfg'.format("userCfg"), default=default,
        help="Config-file (yaml) for {}. (default: {})".format("user-specific configure file", default))

    logger.setLevel(logging.WARN + 10*args.quiet - 10*args.verbose)
    aux_logger.setLevel(logging.WARN + 10*args.quiet - 10*args.verbose)
    if not os.path.exists("phasing"):
        os.mkdir("phasing")

    if not os.path.isdir('/home/users/astar/gis/userrig'):
    # generate the sample configure file according to the output of topMed
    	getSampleConf_cmd = "/mnt/projects/rpd/apps/miniconda3/bin/python  " + LIB_PATH + "/phasing_sample_conf.py --type " + args.seqtype + "  --input ./varCall/out/paste   --sites ./varCall/out/svm  -o ./phasing/sample.yaml";
    else:
	    getSampleConf_cmd = "/home/users/astar/gis/gisshared/rpd/apps/miniconda3/bin/python  " + LIB_PATH + "/phasing_sample_conf.py --type " + args.seqtype + "  --input ./varCall/out/paste   --sites ./varCall/out/svm  -o ./phasing/sample.yaml";    
    print(getSampleConf_cmd);
    os.system(getSampleConf_cmd);    
    sample_cfg = "./phasing/sample.yaml";
    
    if sample_cfg:
        if not os.path.exists(sample_cfg):
            logger.fatal("Config file %s does not exist", sample_cfg)
            sys.exit(1)

    with open(sample_cfg) as fh_cfg:
        yaml_data = yaml.safe_load(fh_cfg)
        unknown_keys = set(yaml_data.keys()) - set(['samples', 'sites','chr'])
        if unknown_keys:
            logger.critical("Found unexpected keys in %s (only 'samples chr'"
                    " and 'sites' allowed): %s", fh_cfg, unknown_keys)
            sys.exit(1)

        chrs, sites, Chr = yaml_data['samples'], yaml_data['sites'], yaml_data['chr']

    # turn arguments into user_data that gets merged into pipeline config
    #
    # generic data first
    user_data = dict()
    user_data['mail_on_completion'] = not args.no_mail
    user_data['sites'] = sites
    user_data['samples'] = chrs
    user_data['chr'] = chrs
    if args.name:
        user_data['analysis_name'] = args.name


    if args.seqtype== "WES":

        pipeline_handler = PipelineHandler(
            PIPELINE_NAME, PIPELINE_BASEDIR, 
            "phasing",user_data,
            Snakefile="Snakefile.split.WES",
            master_q=args.master_q,
            slave_q=args.slave_q,
            params_cfgfile=args.params_cfg,
            modules_cfgfile=args.modules_cfg,
            refs_cfgfile=args.references_cfg,
            cluster_cfgfile=get_cluster_cfgfile(CFG_DIR),
            user_cfgfile=args.userCfg)

    elif args.seqtype== "WGS":
        pipeline_handler = PipelineHandler(
            PIPELINE_NAME, PIPELINE_BASEDIR, 
            "phasing",user_data,
            Snakefile="Snakefile.split.WGS",
            master_q=args.master_q,
            slave_q=args.slave_q,
            params_cfgfile=args.params_cfg,
            modules_cfgfile=args.modules_cfg,
            refs_cfgfile=args.references_cfg,
            cluster_cfgfile=get_cluster_cfgfile(CFG_DIR),
            user_cfgfile=args.userCfg)


    pipeline_handler.setup_env()
    pipeline_handler.submit(args.no_run)

if __name__ == "__main__":
    main()
