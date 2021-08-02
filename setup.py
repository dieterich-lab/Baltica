#!/usr/bin/env python
import pathlib
from setuptools.command.install import install
from setuptools import setup

from baltica.version import __version__, _program

with open('requirements.txt') as f:
    required = [x for x in f.read().splitlines() if not x.startswith("#")]

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

test_required = ['pytest', 'pyfakefs']
# use a custom install https://blog.niteo.co/setuptools-run-custom-code-in-setup-py/
class CustomInstallCommand(install):
    """"""

    def run(self):
        install.run(self)


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
      {program} = baltica.command:cli
      """.format(program=_program),
      install_requires=required,
      test_requires=test_required,
      include_package_data=True,
      scripts=[
          "baltica/annotate_SJ.R",
          "baltica/parse_junctionseq_output.R",
          "baltica/parse_leafcutter_output.R",
          "baltica/parse_majiq_output.R",
          "baltica/assign_AS_type.R",
          "baltica/simplify.R",
          "baltica/parse_rmats_output.R",
          "baltica/baltica_report.Rmd"
      ],
      keywords=['differential splicing', 'bioinformatics', 'rna-seq'],
      zip_safe=False,
      project_urls={
          'Bug Reports': 'https://github.com/dieterich-lab/baltica/issues',
          'Dieterich Lab': 'https://dieterichlab.org',
          'Source': 'https://github.com/dieterich-lab/baltica'
      })
