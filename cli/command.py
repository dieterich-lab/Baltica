#!/usr/bin/env python
"""
Command line interface driver for snakemake workflows
Based on https://github.com/charlesreid1/2019-snakemake-cli
"""
import argparse

import os.path
from pathlib import Path

import snakemake
import sys
import yaml

from . import _program

this_dir = os.path.abspath(os.path.dirname(__file__))
parent_dir = os.path.join(this_dir, '..')
cwd = os.getcwd()

args = tuple(sys.argv[1:])


def main(_args=args):
    parser = argparse.ArgumentParser(prog=_program,
                                     description='Baltica: One stop solution for differential splicing analysis.',
                                     usage='''baltica <workflow> <parameters> [<target>]
Baltica: workflows for alternative splicing analysis.
''')
    parser.add_argument('workflow', required=True)
    parser.add_argument('program', choices=['leafcutter', 'majiq', 'junctionseq'], required=True)
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-f', '--force', action='store_true')
    parser.add_argument('--use-conda', action='store_true')
    parser.add_argument('--use-envmodule', action='store_true')
    _args = parser.parse_args(_args)

    # check if workflow file is readable
    p = Path(__file__) / 'rules' / _args.workflow
    snakefile = p.with_suffix('.smk')

    if not os.path.exists(snakefile):
        msg = 'Error: cannot find Snakefile for the selected program: {} \n'
        sys.stderr.write(msg.format(snakefile))
        sys.exit(-1)
    # TODO check dependencies option

    # check the config file for some mandatory parameters
    parameters = {
        "general": "sample_path workdir samples ref contrasts".split(),
        "leafcutter": "min_samples_per_group min_samples_per_intron fdr".split(),
        "majiq": "threshold GRCh38_90 strandness read_len".split(),
        "junctionseq": "stradness red_len".split()
    }

    with open(snakefile) as fin:
        workflow_info = yaml.safe_load(fin)

    expected_par = set(parameters['general'] + parameters[_args.program])
    missing = expected_par.difference(workflow_info.keys())

    if missing:
        msg = 'Error: the following parameters are missing from the selected workflow file ({}): \n{}'
        sys.stderr.write(msg.format(snakefile, ', '.join(missing)))
        sys.exit(-1)

    target = workflow_info['path']

    print('--------')
    print('details:')
    print('\tsnakefile: {}'.format(snakefile))
    print('\tconfig: {}'.format(snakefile))
    print('\ttarget: {}'.format(target))
    print('--------')

    # run Baltica!
    status = snakemake.snakemake(
        snakefile,
        configfile=_args.params,
        targets=[target],
        printshellcmds=True,
        dryrun=_args.dry_run,
        forceall=_args.force,
        use_conda=_args.use_conda,
        use_envmodule=_args.use_envmodule
    )

    if status:  # translate "success" into shell exit code of 0
        return 0
    return 1


if __name__ == '__main__':
    main()
