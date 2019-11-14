#!/usr/bin/env python

"""
Command line interface driver for snakemake workflows
Based on https://github.com/charlesreid1/2019-snakemake-cli
"""
import argparse
import os.path
import snakemake
import sys
import yaml
from pathlib import Path

from . import _program

this_dir = os.path.abspath(os.path.dirname(__file__))
parent_dir = os.path.join(this_dir, '..')
cwd = os.getcwd()


def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(prog=_program,
                                     description='Baltica: One stop solution for differential splicing analysis.',
                                     usage='''baltica <workflow> <parameters> [<target>]
Baltica: workflows for alternative splicing analysis.
''')
    rules_path = Path('../rules/')
    parser.add_argument('workflowfile', required=True)
    parser.add_argument('paramsfile', required=True)
    parser.add_argument('program', choices=['leafcutter'], required=True)
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-f', '--force', action='store_true')
    args = parser.parse_args(args)

    # first, find the Snakefile
    snakefile = args.program + '.smk'
    if not os.path.exists(snakefile):
        msg = 'Error: cannot find Snakefile for the selected program: {} \n'
        sys.stderr.write(msg.format(args.program))
        sys.exit(-1)

    # next, find the workflow config file
    workflow_file = args.workflowfile
    if not os.path.exists(workflow_file):
        msg = 'Error: cannot find the selected workflow file: {} \n'
        sys.stderr.write(msg.format(workflow_file))
        sys.exit(-1)

    # finally, check the config file for some mandatory parameters
    parameters = {
        "general": "sample_path workdir samples ref contrasts".split(),
        "leafcutter": "min_samples_per_group min_samples_per_intron fdr".split(),
        "majiq": "threshold GRCh38_90 strandness read_len".split()
    }

    with open(workflow_file) as fin:
        workflow_info = yaml.safe_load(fin)

    expected_par = set(parameters['general'] + parameters[args.program])
    missing = expected_par.difference(workflow_info.keys())

    if missing:
        msg = 'Error: the following parameters are missing from the selected workflow file ({}): \n{}'
        sys.stderr.write(msg.format(workflow_file, ', '.join(missing)))
        sys.exit(-1)

    target = workflow_info['workflow_target']

    print('--------')
    print('details')
    print('\tsnakefile: {}'.format(snakefile))
    print('\tconfig: {}'.format(workflow_file))
    print('\ttarget: {}'.format(target))
    print('--------')

    # run Baltica!
    status = snakemake.snakemake(snakefile, configfile=paramsfile,
                                 targets=[target], printshellcmds=True,
                                 dryrun=args.dry_run, forceall=args.force,
                                 config=config)

    if status:  # translate "success" into shell exit code of 0
        return 0
    return 1


if __name__ == '__main__':
    main()
