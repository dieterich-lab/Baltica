#!/usr/bin/env python
"""
Command line interface for Baltica

"""
import os
from pathlib import Path
import sys
import yaml
import logging
import subprocess

import snakemake
import click

_program = "baltica"
__version__ = "1.0.1"

baltica_path = Path(__file__)


# from https://stackoverflow.com/a/56944256/1694714
class CustomFormatter(logging.Formatter):
    """Logging Formatter to add colors and count warning / errors"""

    grey = "\x1b[38;21m"    
    green = "\x1b[32;21m"
    yellow = "\x1b[33;21m"
    red = "\x1b[31;21m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "baltica:\t| %(message)s"

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: green + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

def avaiable_workflows(baltica_path):
    smk_files = baltica_path.parent.glob('*.smk')
    return [x.with_suffix('').name for x in smk_files]

avaiable_workflows_ = avaiable_workflows(baltica_path)

logger = logging.getLogger(__file__)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
logger.addHandler(ch)

# https://click.palletsprojects.com/en/8.0.x/advanced/#forwarding-unknown-options
# unknown options are passed to snakemake
@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument("workflow", type=click.Choice(avaiable_workflows_, case_sensitive=False), default='baltica')
@click.argument("config_file",  type=click.Path(exists=True))
@click.option('-v', '--verbose', is_flag=True, help='Enables verbose mode.')  
@click.argument('snakemake_args', nargs=-1, type=click.UNPROCESSED) 
def baltica(workflow, config_file, verbose, snakemake_args):
    f"""{_program} {__version__}implements workflows for differential junction
     usage methods, and method integration and analysis. Visit
      https://github.com/dieterich-lab/Baltica for more information. 
      
      Runs baltica WORKFLOW with CONFIG_FILE and SNAKEMAKE_ARGS"""
    # TODO add link to baltica docs with important snakemake parameters
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.ERROR)

    # Error if older version of snakemake is installed
    min_snakemake_version = "5.2"
    try:
        snakemake.utils.min_version(min_snakemake_version)
    except snakemake.exceptions.WorkflowError as e:
        logger.error(
            f'{_program} requires Snakemake version >= {min_snakemake_version}:',  
            exc_info=True)
        sys.exit(1)
    # check if workflow file is readable
    snakefile = (baltica_path.parent / workflow).with_suffix(".smk")

    with open(config_file) as fin:
        config = yaml.safe_load(fin)

    target = config["path"]
    
    try:
        os.makedirs(Path(target) / 'logs/')
    except FileExistsError:
        pass

    snakemake_args = list(snakemake_args)
    print(snakemake_args)
    if verbose:
        snakemake_args.extend(['--printshellcmds', '--verbose', '--reason'])
    
    if not any([x in snakemake_args for x in ['--cores', '-c', '--job', '-j']]):
        snakemake_args.append('-j1')

    # here we set bindings three directories needed by singularity
    # the target path, where the output is written to
    # the sample path, which contains the input data
    # the baltica directory, which contains the analysis scripts
    if '--use-singularity' in snakemake_args and "--singularity-args" not in snakemake_args:
        relative_path = os.path.relpath(Path.home(), baltica_path)
        relative_path = Path(relative_path).resolve()
        home_baltica_path = Path.home() / '/'.join(baltica_path.parts[len(relative_path.parts):])
        snakemake_args.extend(['--singularity-args', f'-B {config["path"]},{config["sample_path"]},{home_baltica_path}'])
    
    tmpdir = os.environ['TMPDIR']
    if '--use-singularity' and tmpdir != '/tmp':
        logger.warn("""Current TMPDIR is {tmpdir}, make sure it is writtable 
        by singularity by using --singularity-args or setting TMPDIR to /tmp""")

    # Singularity support 
    try:
        _ = subprocess.run(['singularity', '--version'], stdout=subprocess.DEVNULL)
    except FileNotFoundError as e:
        if '--use-singularity' in snakemake_args:
            logger.critical("Baltica requires singularity, which was not found", exc_info=True)
            sys.exit(1)
    
    logger.info(
        f"""Starting baltica (v{__version__}) analysis with:
    WORKFLOW: {workflow} (from {snakefile})
    CONFIGURATION: {config_file}
    TARGET DIRECTORY: {target}   
    SNAKEMAKE ARGUMENTS: {' '.join(snakemake_args)}
    """)
    cmd = [
        'snakemake',
        '--configfile', config_file,
        '--snakefile', str(snakefile),
        *snakemake_args]
    logger.debug('Start of snakemake logger:')
    result = subprocess.run(cmd)
    return result.check_returncode

if __name__ == "__main__":
    baltica()
