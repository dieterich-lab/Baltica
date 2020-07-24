#!/usr/bin/env python
import os
import pathlib
import sys
from setuptools.command.install import install
from setuptools import setup

from baltica.command import __version__, _program

with open('requirements.txt') as f:
    required = [x for x in f.read().splitlines() if not x.startswith("#")]

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()


def post_install():
    import shutil
    import fileinput

    path = shutil.which("leafcutter_cluster.py")
    with fileinput.input(files=path, inplace=True) as f:
        for i, line in enumerate(f):
            if i != 0:
                print(line, end='')
            else:
                print('#!/usr/bin/env python2')


# use a custom install https://blog.niteo.co/setuptools-run-custom-code-in-setup-py/
class CustomInstallCommand(install):
    """Customized setuptools install command - prints a friendly greeting."""
    def run(self):
        install.run(self)
        print('Running post install task')
        post_install()



setup(name=_program,
      cmdclass={'install': CustomInstallCommand},
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
          "scripts/leafcutter_ds_pair.R",
          "scripts/parse_junctionseq_output.R",
          "scripts/parse_leafcutter_output.R",
          "scripts/parse_majiq_output.R",
          "scripts/sam2bed.pl",
          "scripts/utils.R",
          "scripts/assign_AS_type.R",
          "scripts/simplify.R"

          ],
      keywords=['differential splicing', 'bioinformatics', 'rna-seq'],
      zip_safe=False,
      project_urls={
          'Bug Reports': 'https://github.com/dieterich-lab/baltica/issues',
          'Dieterich Lab': 'https://dieterichlab.org',
          'Source': 'https://github.com/dieterich-lab/baltica'
      })
