"""library functions for pipelines
"""

#--- standard library imports
#
import os
import sys
import subprocess
import logging
import shutil
import smtplib
#from email.mime.text import MIMEText
from getpass import getuser
#import socket
import time
from datetime import datetime, timedelta
import calendar
import json
import requests
from pathlib import Path

#--- third-party imports
#
import yaml

#--- project specific imports
#
#from services import SMTP_SERVER
#from services import rest_services




# only dump() and following do not automatically create aliases
yaml.Dumper.ignore_aliases = lambda *args: True


# global logger
logger = logging.getLogger(__name__)
# handler = logging.StreamHandler()
# handler.setFormatter(logging.Formatter(
#     '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
# logger.addHandler(handler)

# dir relative to Snakefile where configs are to be found
CFG_DIR = "cfg"

PIPELINE_ROOTDIR = os.path.join(os.path.dirname(__file__), "..")
assert os.path.exists(os.path.join(PIPELINE_ROOTDIR, "VERSION"))

class PipelineHandler(object):
    """FIXME:add-doc

    - FIXME needs cleaning up!
    - FIXME check access of global vars
    """

    # output
    PIPELINE_CFGFILE = "conf.yaml"

    RC_DIR = "rc"

    RC_FILES = {
        'DK_INIT' : os.path.join(RC_DIR, 'dk_init.rc'),# used to load dotkit
        'SNAKEMAKE_INIT' : os.path.join(RC_DIR, 'snakemake_init.rc'),# used to load snakemake
        'SNAKEMAKE_ENV' : os.path.join(RC_DIR, 'snakemake_env.rc'),# used as bash prefix within snakemakejobs
    }

    LOG_DIR_REL = "logs"
    # MASTERLOG = os.path.join(LOG_DIR_REL, "snakemake.log")
    SUBMISSIONLOG = os.path.join(LOG_DIR_REL, "submission.log")

    # master max walltime in hours
    # note, this includes waiting for jobs in q
    MASTER_WALLTIME_H = 120

    def __init__(self,
                 pipeline_name,
                 pipeline_subdir,
                 outdir,
                 user_data,
                 pipeline_rootdir=PIPELINE_ROOTDIR,
                 logger_cmd="true", # bash: not doing anything by default
                 site=None,
                 Snakefile=None,
                 master_q=None,
                 slave_q=None,
                 master_walltime_h=MASTER_WALLTIME_H,
                 params_cfgfile=None,
                 modules_cfgfile=None,
                 refs_cfgfile=None,
                 user_cfgfile=None,
                 cluster_cfgfile=None,
                 local_mode=False):
        """FIXME:add-doc

        pipeline_subdir: where default configs can be found, i.e pipeline subdir
        """

        self.outdir = outdir
        self.pipeline_name = pipeline_name
        self.pipeline_version = get_pipeline_version()# external function
        self.pipeline_subdir = pipeline_subdir
        self.user_data = user_data

        self.log_dir_rel = self.LOG_DIR_REL

        # make the masterlog file name more informative
        self.masterlog = pipeline_name+".master.log"
        self.submissionlog = self.SUBMISSIONLOG
        self.local_mode=local_mode


        #if params_cfgfile:
        #    assert os.path.exists(params_cfgfile)
        #self.params_cfgfile = params_cfgfile

        #if modules_cfgfile:
        #    assert os.path.exists(modules_cfgfile)
        #self.modules_cfgfile = modules_cfgfile

        #if refs_cfgfile:
        #    assert os.path.exists(refs_cfgfile)
        #self.refs_cfgfile = refs_cfgfile

        if user_cfgfile:
            assert os.path.exists(user_cfgfile)
        self.user_cfgfile = user_cfgfile

        if cluster_cfgfile:
            assert os.path.exists(cluster_cfgfile)
        self.cluster_cfgfile = cluster_cfgfile

        self.pipeline_cfgfile_out = os.path.join(
            self.outdir, self.PIPELINE_CFGFILE)

        # RCs
        self.dk_init_file = os.path.join(
            self.outdir, self.RC_FILES['DK_INIT'])
        self.snakemake_init_file = os.path.join(
            self.outdir, self.RC_FILES['SNAKEMAKE_INIT'])
        self.snakemake_env_file = os.path.join(
            self.outdir, self.RC_FILES['SNAKEMAKE_ENV'])

        self.logger_cmd = logger_cmd
        self.master_q = master_q
        self.slave_q = slave_q
        self.master_walltime_h = master_walltime_h
        self.snakefile_abs = os.path.abspath(
            os.path.join(pipeline_subdir, Snakefile))
        assert os.path.exists(self.snakefile_abs)

        # cluster configs
        if self.cluster_cfgfile:
            self.cluster_cfgfile_out = os.path.join(outdir, "cluster.yaml")
        # else: local

        # run template
        self.run_template = os.path.join(
            pipeline_rootdir, "lib", "run.template.sh")
        self.run_out = os.path.join(outdir, "run.sh")
        assert os.path.exists(self.run_template)

        log_path = os.path.abspath(os.path.join(self.outdir, self.masterlog))
        self.elm_data = {'pipeline_name': self.pipeline_name,
                         'pipeline_version': self.pipeline_version,
                         'instance_id': 'SET_ON_EXEC',# dummy
                         'submitter': 'SET_ON_EXEC',# dummy
                         'log_path': log_path}


    @staticmethod
    def write_dk_init(rc_file, overwrite=True):
        """write dotkit init rc
        """
        if not overwrite:
            assert not os.path.exists(rc_file), rc_file
        with open(rc_file, 'w') as fh:
            fh.write("eval `{}`;\n".format(' '.join(get_init_call())))


    @staticmethod
    def write_snakemake_init(rc_file, overwrite=True):
        """write snakemake init rc (loads miniconda and, activate source')
        """
        if not overwrite:
            assert not os.path.exists(rc_file), rc_file
        with open(rc_file, 'w') as fh:
            fh.write("# initialize snakemake. requires pre-initialized dotkit\n")
            #fh.write("reuse -q miniconda-3\n")
            #fh.write("source activate snakemake-3.5.5-g9752cd7-catch-logger-cleanup\n")
            #fh.write("source activate snakemake-4.8\n")


    def write_snakemake_env(self, overwrite=True):
        """creates rc file for use as 'bash prefix', which also loads modules defined in cfgfile
        """

        if not overwrite:
            assert not os.path.exists(self.snakemake_env_file), self.snakemake_env_file

        with open(self.snakemake_env_file, 'w') as fh_rc:
            fh_rc.write("# used as bash prefix within snakemake\n\n")
            fh_rc.write("# init dotkit\n")
            fh_rc.write("source {}\n\n".format(os.path.relpath(self.dk_init_file, self.outdir)))

            fh_rc.write("# load modules\n")
            #with open(self.pipeline_cfgfile_out) as fh_cfg:
            #    yaml_data = yaml.safe_load(fh_cfg)
            #    assert "modules" in yaml_data
            #    for k, v in yaml_data["modules"].items():
            #        #fh_rc.write("reuse -q {}\n".format("{}-{}".format(k, v)))

            fh_rc.write("\n")
            fh_rc.write("# unofficial bash strict has to come last\n")
            fh_rc.write("set -euo pipefail;\n")


    def write_cluster_config(self):
        """writes site dependend cluster config
        """
        shutil.copyfile(self.cluster_cfgfile, self.cluster_cfgfile_out)


    def write_run_template(self):
        """FIXME:add-doc
        """

        # have to specify the name of working directory so that the main script
        # can run properly on Torque
        d = {'SNAKEFILE': self.snakefile_abs,
             'LOGDIR': self.log_dir_rel,
             'MASTERLOG': self.log_dir_rel+"/"+self.masterlog,
             'PIPELINE_NAME': self.pipeline_name,
             'MASTER_WALLTIME_H': self.master_walltime_h,
             'DEFAULT_SLAVE_Q': self.slave_q if self.slave_q else "",
             'LOGGER_CMD': self.logger_cmd,
	         'PIPELINE_ROOT': PIPELINE_ROOTDIR,
             'WORKING_DIR': os.getcwd()+"/"+self.outdir}

        with open(self.run_template) as fh:
            templ = fh.read()
        with open(self.run_out, 'w') as fh:
            fh.write(templ.format(**d))


    def read_cfgfiles(self):
        """parse default config and replace all RPD env vars
        """

        merged_cfg = dict()
        #rpd_vars = get_rpd_vars()

        for cfgkey, cfgfile in [
                                ('users', self.user_cfgfile),
                                ]:
            if not cfgfile:
                continue
            with open(cfgfile) as fh:
                cfg = dict(yaml.safe_load(fh))
            # to replace rpd vars the trick is to traverse
            # dictionary fully and replace all instances
            dump = json.dumps(cfg)
            #for k, v in rpd_vars.items():
            #    dump = dump.replace("${}".format(k), v)
            #cfg = dict(json.loads(dump))
            if cfgkey == 'global':
                merged_cfg.update(cfg)
            else:
                assert cfgkey not in merged_cfg
                merged_cfg[cfgkey] = cfg

        # determine num_chroms needed by some pipelines
        # FIXME ugly because sometimes not needed
        #if merged_cfg.get('references'):
        #    reffa = merged_cfg['references']['genome']
        #    if reffa:
        #        assert 'num_chroms' not in merged_cfg['references']
        #        merged_cfg['references']['num_chroms'] = len(list(
        #            chroms_and_lens_from_from_fasta(reffa)))

        return merged_cfg



    def write_merged_cfg(self, force_overwrite=True):
        """writes config file for use in snakemake becaused on default config
        """

        config = self.read_cfgfiles()
        config.update(self.user_data)
        assert 'ELM' not in config
        config['ELM'] = self.elm_data

        print(self.pipeline_cfgfile_out)
        if not force_overwrite:
            assert not os.path.exists(self.pipeline_cfgfile_out)
        with open(self.pipeline_cfgfile_out, 'w') as fh:
            # default_flow_style=None(default)|True(least readable)|False(most readable)
            yaml.dump(config, fh, default_flow_style=False)


    def setup_env(self):
        """FIXME:add-doc
        """

        logger.info("Creating run environment in %s", self.outdir)
        # create log dir recursively so that parent is created as well
        if not os.path.exists(os.path.join(self.outdir, self.log_dir_rel)):
            os.makedirs(os.path.join(self.outdir, self.log_dir_rel))
            os.makedirs(os.path.join(self.outdir, self.RC_DIR))

        if not self.local_mode:
            self.write_cluster_config()
        self.write_merged_cfg()
        self.write_snakemake_env()
        self.write_dk_init(self.dk_init_file)
        self.write_snakemake_init(self.snakemake_init_file)
        self.write_run_template()



    def submit(self, no_run=False):
        """FIXME:add-doc
        """

        if self.master_q:
            master_q_arg = "-q {}".format(self.master_q)
        else:
            master_q_arg = ""
        if self.local_mode:
            logger.warning("Please note that script is run in 'local' mode"
                           " (which is mainly for debugging)")
            cmd = "cd {} && bash {} {} >> {}".format(
                os.path.dirname(self.run_out), master_q_arg,
                os.path.basename(self.run_out), self.submissionlog)
        else:
            cmd = "cd {} && qsub {} {} >> {}".format(
                os.path.dirname(self.run_out), master_q_arg,
                os.path.basename(self.run_out), self.submissionlog)

        if no_run:
            logger.warning("Skipping pipeline run on request. Once ready, use: %s", cmd)
            logger.warning("Once ready submit with: %s", cmd)
        else:
            logger.info("Starting pipeline: %s", cmd)
            #os.chdir(os.path.dirname(run_out))
            _ = subprocess.check_output(cmd, shell=True)
            submission_log_abs = os.path.abspath(os.path.join(self.outdir, self.submissionlog))
            master_log_abs = os.path.abspath(os.path.join(self.outdir, self.masterlog))
            logger.debug("For submission details see %s", submission_log_abs)
            logger.info("The (master) logfile is %s", master_log_abs)



def get_pipeline_version():
    """determine pipeline version as defined by updir file
    """
    version_file = os.path.abspath(os.path.join(PIPELINE_ROOTDIR, "VERSION"))
    with open(version_file) as fh:
        version = fh.readline().strip()
    cwd = os.getcwd()
    os.chdir(PIPELINE_ROOTDIR)
    if os.path.exists(".git"):
        commit = None
        cmd = ['git', 'describe', '--always', '--dirty']
        try:
            res = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            commit = res.decode().strip()
        except (subprocess.CalledProcessError, OSError) as _:
            pass
        if commit:
            version = "{} commit {}".format(version, commit)
    os.chdir(cwd)
    return version


def is_devel_version():
    """checks whether this is a developers version of production
    """
    check_file = os.path.abspath(os.path.join(PIPELINE_ROOTDIR, "DEVELOPERS_VERSION"))
    #logger.debug("check_file = {}".format(check_file))
    return os.path.exists(check_file)


def get_cluster_cfgfile(cfg_dir):
    """returns None for local runs
    """
    cfg = os.path.join(cfg_dir, "cluster.yaml")
    assert os.path.exists(cfg), ("Missing file {}".format(cfg))
    return cfg


def get_init_call():
    """return dotkit init call
    """
    path_to_init=os.path.join(PIPELINE_ROOTDIR, "init")
    if not os.path.exists(path_to_init):
        raise Exception("Initialization file "+path_to_init+" not found!")

    return path_to_init


def get_rpd_vars():
    """Read RPD variables set by calling and parsing output from init
    """

    cmd = get_init_call()
    try:
        res = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        logger.fatal("Couldn't call init as '%s'", ' '.join(cmd))
        raise

    rpd_vars = dict()
    for line in res.decode().splitlines():
        if line.startswith('export '):
            line = line.replace("export ", "")
            line = ''.join([c for c in line if c not in '";\''])
            #logger.debug("line = {}".format(line))
            k, v = line.split('=')
            rpd_vars[k.strip()] = v.strip()
    return rpd_vars



def generate_timestamp():
    """generate ISO8601 timestamp incl microsends, but with colons
    replaced to avoid problems if used as file name
    """
    return datetime.isoformat(datetime.now()).replace(":", "-")


def timestamp_from_string(analysis_id):
    """
    converts output of generate_timestamp(), e.g. 2016-05-09T16-43-32.080740 back to timestamp
    """
    dt = datetime.strptime(analysis_id, '%Y-%m-%dT%H-%M-%S.%f')
    return dt


def isoformat_to_epoch_time(ts):
    """
    Converts ISO8601 format (analysis_id) into epoch time
    """
    dt = datetime.strptime(ts[:-7], '%Y-%m-%dT%H-%M-%S.%f')-\
         timedelta(hours=int(ts[-5:-3]),
                   minutes=int(ts[-2:]))*int(ts[-6:-5]+'1')
    epoch_time = calendar.timegm(dt.timetuple()) + dt.microsecond/1000000.0
    return epoch_time


def get_machine_run_flowcell_id(runid_and_flowcellid):
    """return machine-id, run-id and flowcell-id from full string.
    Expected string format is machine-runid_flowcellid
    """
    # be lenient and allow full path
    runid_and_flowcellid = runid_and_flowcellid.rstrip("/").split('/')[-1]

    runid, flowcellid = runid_and_flowcellid.split("_")
    machineid = runid.split("-")[0]
    return machineid, runid, flowcellid


def ref_is_indexed(ref, prog="bwa"):
    """checks whether a reference was already indexed by given program"""

    if prog == "bwa":
        return all([os.path.exists(ref + ext)
                    for ext in [".pac", ".ann", ".amb", ".sa"]])
    elif prog == "samtools":
        return os.path.exists(ref + ".fai")
    else:
        raise ValueError


def generate_window(days=7):
    """returns tuple representing epoch window (int:present, int:past)"""
    date_time = time.strftime('%Y-%m-%d %H:%M:%S')
    pattern = '%Y-%m-%d %H:%M:%S'
    epoch_present = int(time.mktime(time.strptime(date_time, pattern)))*1000
    d = datetime.now() - timedelta(days=days)
    f = d.strftime("%Y-%m-%d %H:%m:%S")
    epoch_back = int(time.mktime(time.strptime(f, pattern)))*1000
    return (epoch_present, epoch_back)


def parse_regions_from_bed(bed):
    """yields regions from bed as three tuple
    """

    with open(bed) as fh:
        for line in fh:
            if line.startswith('#') or not len(line.strip()) or line.startswith('track '):
                continue
            chrom, start, end = line.split()[:3]
            start, end = int(start), int(end)
            yield (chrom, start, end)


def chroms_and_lens_from_from_fasta(fasta):
    """return sequence and their length as two tuple. derived from fai
    """

    fai = fasta + ".fai"
    assert os.path.exists(fai), ("{} not indexed".format(fasta))
    with open(fai) as fh:
        for line in fh:
            (s, l) = line.split()[:2]
            l = int(l)
            yield (s, l)
