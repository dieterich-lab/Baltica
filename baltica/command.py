#!/usr/bin/env python
"""
Command line interface driver for Baltica

Modified from https://github.com/charlesreid1/2019-snakemake-cli
"""
import os
import argparse
from pathlib import Path
import sys
import yaml

import snakemake

try:
    import baltica
    baltica_installed = True
except ImportError:
    baltica_installed = False

_program = "baltica"
__version__ = "1.0.1"
p = Path(__file__)

desc = f"{_program} implements workflows for differential junction usage and consequence analysis. Visit " \
       f"https://github.com/dieterich-lab/Baltica for more information. "


def main():
    parser = argparse.ArgumentParser(prog=_program,
                                     description=desc,
                                     usage=f"{_program} <workflow> <config> <options>")
    parser.add_argument(
        "workflow",
        choices=["leafcutter", "majiq", "junctionseq", "qc", "stringtie", "analysis"],
        help="Workflow to be run"
    )
    parser.add_argument(
        "config",
        type=Path,
        help="Configuration file for the workflow."
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Execute a test run and output commands to be run."
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force the workflow to run despite the presence of an output."
    )
    parser.add_argument(
        "--use-conda",
        action="store_true",
        help="Use conda to install the requirements."
    )
    # TODO add slurm profile
    # parser.add_argument(
    #     "--profile",
    #     default="slurm",
    #     help="Snakemake cluster profile used to run the workflows"
    # ),
    parser.add_argument(
        "--use-envmodule",
        action="store_true",
        help="Use environment modules for requirements."
    )
    parser.add_argument(
        "--unlock",
        action="store_true",
        help="Use directory from previously running DAG that was interrupted."
    ),
    parser.add_argument(
        "--cluster",
        metavar='str',
        default="sbatch --parsable --mem {cluster.mem} --out {cluster.out} --error {cluster.out} -c {cluster.cpu} --partition {cluster.partition}",
        help='Snakemake cluster parameter (default: %(default)s)'
    ),
    parser.add_argument(
        "--cluster-config",
        metavar='str',
        default=str(p.parent / 'cluster.yml'),
        help='Snakemake cluster configuration (default: %(default)s)'
    ),
    parser.add_argument(
        "--nodes",
        type=int,
        default=10
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__
    )

    # Error if older version of snakemake is installed
    min_snakemake_version = "5.2"
    try:
        snakemake.utils.min_version(min_snakemake_version)
    except snakemake.exceptions.WorkflowError as e:
        print(f'{_program} requires Snakemake version >= {min_snakemake_version}:', file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(1)

    argv = parser.parse_args(tuple(sys.argv[1:]))
    # check if workflow file is readable
    snakefile = (p.parent / argv.workflow).with_suffix(".smk")

    with open(argv.config) as fin:
        workflow_info = yaml.safe_load(fin)

    target = workflow_info["path"]

    if not snakefile.exists:
        print(f"Error: cannot find Snakefile for the {argv.workflow} workflow: {snakefile} \n", file=sys.stderr)
        sys.exit(1)

    if not argv.config.exists:
        print(f"Error: cannot find configuration file {argv.config} \n", file=sys.stderr)
        sys.exit(1)

    # TODO handling cluster profile profile_config = Path(snakemake.get_profile_file(argv.profile, 'config.yaml',
    #  return_default=True)) SBATCH_DEFAULTS = """'sbatch -p {cluster.partition} --mem {cluster.mem} --out {
    #  cluster.out} --error {cluster.out} -c {cluster.cpu}'""" CLUSTER_CONFIG = "cluster.yaml"

    print("--------")
    print("details:")
    print("\tsnakefile: {}".format(snakefile))
    print("\tconfiguration file: {}".format(argv.config))
    print("\twork directory: {}".format(target))
    print("--------")

    try:
        os.makedirs(Path(target) / 'logs/')
    except FileExistsError:
        pass

    # run Baltica workflow
    success = snakemake.snakemake(
        snakefile,
        configfiles=[argv.config],
        workdir=str(target),
        printshellcmds=True,
        dryrun=argv.dry_run,
        forceall=argv.force,
        use_conda=argv.use_conda,
        use_env_modules=argv.use_envmodule,
        cluster=argv.cluster,
        cluster_config=argv.cluster_config,
        cluster_sync=None,
        unlock=argv.unlock,
        nodes=argv.nodes
    )

    if success:
        return 0
    return 1


if __name__ == "__main__":
    main()
