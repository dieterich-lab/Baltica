#!/usr/bin/env python

from setuptools import setup
from baltica.command import __version__, _program
import pathlib

with open('requirements.txt') as f:
    required = [x for x in f.read().splitlines() if not x.startswith("#")]

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

setup(name=_program,
      version=__version__,
      packages=['baltica'],
      description='Workflows for differential junction usage with Baltica',
      long_description=README,
      url='https://github.com/dieterich-lab/Baltica',
      author='@tbrittoborges',
      author_email='thiago.brittoborges@uni-heidelberg.de',
      license='MIT',
      classifiers=[
          "License :: OSI Approved :: MIT License",
          "Programming Language :: Python :: 3",
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          "Intended Audience :: Science/Research"],
      entry_points="""
      [console_scripts]
      {program} = baltica.command:main
      """.format(program=_program),
      install_requires=required,
      include_package_data=True,
      scripts=[
          "scripts/annotate_SJ.R",
          "scripts/bam2junc.sh",
          "scripts/bed2junc.pl",
          "scripts/ds_plots.R",
          "scripts/filter_cs.py",
          "scripts/gtf2gff3.pl",
          "scripts/gtf_to_exons.R",
          "scripts/junctionSeq.R",
          "scripts/leafcutter_cluster.py",
          "scripts/leafcutter_ds.R",
          "scripts/parse_junctionseq_output.R",
          "scripts/parse_leafcutter_output.R",
          "scripts/parse_majiq_output.R",
          "scripts/sam2bed.pl",
          "scripts/utils.R"],
      keywords=['differential splicing', 'bioinformatics', 'rna-seq'],
      zip_safe=False,
      project_urls={
          'Bug Reports': 'https://github.com/dieterich-lab/baltica/issues',
          'Dieterich Lab': 'https://dieterichlab.org',
          'Source': 'https://github.com/dieterich-lab/baltica'
      })
