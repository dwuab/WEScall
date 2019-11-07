#!/usr/bin/env python3
"""Creates a config file describing your samples that can be used as
input for phasing pipeline (-i)

"""

#--- standard library imports
#
import sys
import os
import re
import argparse
import logging
import glob
import csv

#--- third-party imports
#
import yaml

#--- project specisfic imports
#
# add lib dir for this pipeline installation to PYTHONPATH
LIB_PATH = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
if LIB_PATH not in sys.path:
    sys.path.insert(0, LIB_PATH)
#from readunits import ReadUnit
#from readunits import create_rg_id_from_ru
#from readunits import key_for_readunit

__author__ = "Jinzhuang Dou"
__email__ = "douj@gis.a-star.edu.sg"
__copyright__ = "2016 Genome Institute of Singapore"
__license__ = "The MIT License (MIT)"


# only dump() and following do not automatically create aliases
yaml.Dumper.ignore_aliases = lambda *args: True


# global logger
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)


def main():
    """main function
    """

    parser = argparse.ArgumentParser(description=__doc__)

    # generic args
    parser.add_argument('-i', "--input", required=True,
                        help="the directory where your samples are")
    parser.add_argument('-o', "--yaml", required=True,
                        help="Output config (yaml) file")
    parser.add_argument('-f', '--force-overwrite', action='store_true',
                        help="Force overwriting of existing file")
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="Increase verbosity")
    parser.add_argument('-q', '--quiet', action='count', default=0,
                        help="Decrease verbosity")

    args = parser.parse_args()

    logger.setLevel(logging.WARN + 10*args.quiet - 10*args.verbose)

    if not os.path.exists(args.input):
        logger.fatal("Input file %s does not exist", args.input)
        sys.exit(1)
    if os.path.exists(args.yaml) and not args.force_overwrite:
        logger.fatal("Cowardly refusing to overwrite existing file %s", args.yaml)
        sys.exit(1)

    with open(args.input, 'r') as fh:   
        data = yaml.load(fh)
    CHRS = str(data['users']['chrs']).split('-')


    chr_split=dict()
    # Fix it to use the abs path.
    # detect the finished chromosome according to *.OK. 
    for chrFlag in CHRS:  
        for fileName in glob.glob(os.path.join("./phasing/",chrFlag,'*.Split.*.vcf.gz')):  

            chr_split.setdefault(chrFlag, []).append(re.findall(r'\d+', fileName)[-1])

    with open(args.yaml, 'w') as fh:
        yaml.dump(dict(chr_split=chr_split), fh, default_flow_style=False)


if __name__ == "__main__":
    main()
