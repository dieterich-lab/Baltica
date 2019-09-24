from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = [x for x in f.read().splitlines() if not x.startswith("#")]

# Note: the _program variable is set in __init__.py.
# it determines the name of the package/final command line tool.
from cli import __version__, _program

setup(name=_program,
      version=__version__,
      packages=['cli'],
      description='Baltica command line interface',
      url='https://github.com/dieterich-lab/Baltica',
      author='@tbrittoborges',
      author_email='thiago.brittoborges@uni-heidelberg.de',
      license='MIT',
      entry_points="""
      [console_scripts]
      {program} = cli.command:main
      """.format(program=_program),
      install_requires=required,
      include_package_data=True,
      keywords=[],
      zip_safe=False)
