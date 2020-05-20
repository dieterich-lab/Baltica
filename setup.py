#!/usr/bin/env python

from setuptools import setup
from cli import __version__, _program
import pathlib

with open('requirements.txt') as f:
    required = [x for x in f.read().splitlines() if not x.startswith("#")]

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

setup(name=_program,
      version=__version__,
      packages=['cli'],
      description='Workflows for differential splicing with Baltica ',
      longdescription=README,
      url='https://github.com/dieterich-lab/Baltica',
      author='@tbrittoborges',
      author_email='thiago.brittoborges@uni-heidelberg.de',
      license='MIT',
      classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research"],
      entry_points="""
      [console_scripts]
      {program} = cli.command:main
      """.format(program=_program),
      install_requires=required,
      include_package_data=True,
      keywords=[],
      zip_safe=False)
