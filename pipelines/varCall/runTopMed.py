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
import shutil

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

__author__ = "Jinzhuang Dou and others"
__email__ = "douj@gis.a-star.edu.sg"
__copyright__ = "2016 Genome Institute of Singapore"
__license__ = "The MIT License (MIT)"


# only dump() and following do not automatically create aliases
yaml.Dumper.ignore_aliases = lambda *args: True


PIPELINE_BASEDIR = os.path.dirname(sys.argv[0])
CFG_DIR = os.path.join(PIPELINE_BASEDIR, "cfg")

# same as folder name. also used for cluster job names
PIPELINE_NAME = "runTopMed"

# global logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
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
    parser.add_argument('--name',
                        help="Give this analysis run a name (used in email and report)")
    parser.add_argument('--no-mail', action='store_true',
                        help="Don't send mail on completion")
    parser.add_argument('-w', '--slave-q', default=None,
                        help="Queue to use for slave jobs (default: {})".format(None))
    parser.add_argument('-m', '--master-q', default=None,
                        help="Queue to use for master job (default: {})".format(None))
    parser.add_argument('-n', '--no-run', action='store_true')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="Increase verbosity")
    parser.add_argument('-q', '--quiet', action='count', default=0,
                        help="Decrease verbosity")
    parser.add_argument('-i', '--index', required=True, default=0,
                        help="The list of BAM/CRAM files of study samples")
    parser.add_argument('-c', '--userCfg', required=True, default=0,
                        help="User specific configure file")
    parser.add_argument('-t', "--seqtype", required=True,
                        choices=['WGS', 'WES'],
                        help="Sequencing type")
    parser.add_argument('-l', "--local",action='store_true',
        help="Whether to run in local mode")
    cfg_group = parser.add_argument_group('Configuration files (advanced)')
    for name, descr in [("params", "parameters"),
                        ("references", "reference sequences"),
                        ("modules", "modules")]:
        default = os.path.abspath(os.path.join(CFG_DIR, "{}.yaml".format(name)))
        cfg_group.add_argument('--{}-cfg'.format(name),
                               default=default,
                               help="Config-file (yaml) for {}. (default: {})".format(descr, default))

    # Add user-specific configure information
    # 
    # 

    args = parser.parse_args()
    default = os.path.abspath(args.userCfg)

    cfg_group.add_argument('--{}-cfg'.format("userCfg"), default=default,
        help="Config-file (yaml) for {}. (default: {})".format("user-specific configure file", default))


    logger.setLevel(logging.WARN + 10*args.quiet - 10*args.verbose)
    aux_logger.setLevel(logging.WARN + 10*args.quiet - 10*args.verbose)



    # turn arguments into user_data that gets merged into pipeline config
    #
    # generic data first
    user_data = dict()
    user_data['mail_on_completion'] = not args.no_mail
    if args.name:
        user_data['analysis_name'] = args.name

    if not args.local:
        cluster_cfgfile=get_cluster_cfgfile(CFG_DIR)

    pipeline_handler = PipelineHandler(
        PIPELINE_NAME, PIPELINE_BASEDIR, 
        "./varCall",user_data,
        Snakefile="Snakefile."+args.seqtype,
        master_q=args.master_q,
        slave_q=args.slave_q,
        params_cfgfile=args.params_cfg,
        modules_cfgfile=args.modules_cfg,
        refs_cfgfile=args.references_cfg,
        cluster_cfgfile=cluster_cfgfile,
        user_cfgfile=args.userCfg,
        local_mode=args.local)

    ####### 
    os.system("mkdir -p ./varCall/data")
    shutil.copy2(args.index,"./varCall/data/test.index")

    # automatically generate the pedigree file for the user
    with open(args.index) as f_in, open("./varCall/data/test.ped","w") as f_out:
        for line in f_in:
            record = line.strip().split("\t")
            f_out.write("{smp}\t{smp}\t{smp}\t0\t0\n".format(smp=record[0]))

    pipeline_handler.setup_env()
    pipeline_handler.submit(args.no_run)


if __name__ == "__main__":
    main()
