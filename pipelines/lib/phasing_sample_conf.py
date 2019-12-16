#!/usr/bin/env python3
"""Creates a config file describing your samples that can be used as
input for phasing pipeline (-i)

"""

#--- standard library imports
#
import sys
import os
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
    parser.add_argument('-s', "--sites", required=True,
                        help="the directory where the site files are")
    parser.add_argument('-t', "--type", required=True,
                        help="the data type (WES|WGS)")
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

    samples = dict()
    sites = dict()
    chrs=dict()
    # Fix it to use the abs path.
    # detect the finished chromosome according to *.OK. 
    for chrFlag in glob.glob(os.path.join(args.input, '*.OK')):  
        chrFlag=chrFlag.strip('.OK').split('/')[-1]
        for fileName in glob.glob(os.path.join(args.input, chrFlag,'*.bcf')):  
            samples.setdefault(chrFlag, []).append(os.path.abspath(fileName))

    #  fix me
    #  detect the site files after SVM filter. 
    #  now only handle the format "22_1_51304566_milk_transfer.uniq.sites.vcf.gz"
    #  if no site files available, all sites will be assigned to the flag PASS.

    if args.type == 'WGS':
        
        for fullName in glob.glob(os.path.join(args.sites, '*milk_transfer.sites.vcf.gz')):
            fileName=fullName.split('/')[-1]  
            chrFlag=fileName.split('_')[0].split('/')[-1]
            sites.setdefault(chrFlag, []).append(os.path.abspath(fullName))
            chrs.setdefault(chrFlag, []).append("")
    else:
        fullName = os.path.join(args.sites, '0_1_0_milk_svm.sites.vcf.gz');
        for  chrFlag in range(1, 23):
             sites.setdefault(chrFlag, []).append(os.path.abspath(fullName))
             chrs.setdefault(chrFlag, []).append("")



    with open(args.yaml, 'w') as fh:
        yaml.dump(dict(samples=samples), fh, default_flow_style=False)
        yaml.dump(dict(sites=sites), fh, default_flow_style=False)
        yaml.dump(dict(chr=chrs), fh, default_flow_style=False)

if __name__ == "__main__":
    main()
