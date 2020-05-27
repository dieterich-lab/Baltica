#!/usr/bin/env python
"""
Command line interface driver for Baltica

Modified from https://github.com/charlesreid1/2019-snakemake-cli
"""
import argparse
import os.path
from pathlib import Path
import sys
import yaml

import snakemake

# from . import __version__, _program
_program = 'B'
__version__ = '0.0'


def main(_args):
    parser = argparse.ArgumentParser(prog=_program,
                                     description='Baltica: workflows for differential junction usage and '
                                                 'consequence analysis.',
                                     usage='''Baltica <workflow> <config> 
Baltica: workflows for alternative splicing analysis.
''')
    parser.add_argument(
        'workflow', choices=[
            'leafcutter', 'majiq', 'junctionseq', 'qc', 'stringtie', 'analysis'],
        help='Workflow to be run')
    parser.add_argument('config', help='Configuration file for the workflow', type=Path)
    parser.add_argument('-n', '--dry-run', action='store_true', help='')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Force the workflow to run despite the presence of an output')
    parser.add_argument('--use-conda', action='store_true',
                        help='Use conda to install the requirements')
    parser.add_argument('--use-envmodule', action='store_true', help='Use environment modules for requirements')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    _args = parser.parse_args(_args)
    print(_args.config)
    # check if workflow file is readable
    p = Path(__file__).parent.parent / 'rules' / _args.workflow
    snakefile = p.with_suffix('.smk')

    with open(_args.config) as fin:
        workflow_info = yaml.safe_load(fin)

    target = workflow_info['path']

    if not snakefile.exists:
        msg = f'Error: cannot find Snakefile for the {_args.workflow} workflow: {snakefile} \n'
        sys.stderr.write(msg)
        sys.exit(-1)

    if not _args.config.exists:
        msg = f'Error: cannot find configuration file {_args.config} \n'
        sys.stderr.write(msg)
        sys.exit(-1)

    print('--------')
    print('details:')
    print('\tsnakefile: {}'.format(snakefile))
    print('\tconfiguration file: {}'.format(_args.config))
    print('\twork directory: {}'.format(target))
    print('--------')

    # run Baltica workflow
    status = snakemake.snakemake(
        snakefile,
        configfiles=[_args.config],
        workdir=str(target),
        printshellcmds=True,
        dryrun=_args.dry_run,
        forceall=_args.force,
        use_conda=_args.use_conda,
        use_env_modules=_args.use_envmodule
    )

    if status:
        return 0
    return 1


if __name__ == '__main__':
    main(_args=tuple(sys.argv[1:]))
