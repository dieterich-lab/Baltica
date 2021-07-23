from pathlib import Path
from click.testing import CliRunner
from baltica.command import cli, avaiable_workflows
from unittest.mock import MagicMock, patch

def test_baltica_cli_help():
    runner = CliRunner()
    result = runner.invoke(cli, '--help')
    assert result.exit_code == 0

def test_baltica_cli_version():
    runner = CliRunner()
    result = runner.invoke(cli, '--version')
    assert result.exit_code == 0

def test_baltica_cli_version_version():
    runner = CliRunner()
    result = runner.invoke(cli, '--version --version')
    assert result.exit_code == 0

def test_baltica_cli_no_config():
    runner = CliRunner()
    result = runner.invoke(cli, 'rmats')
    assert result.exit_code == 2
    assert result.output.endswith("Error: Missing argument 'CONFIG_FILE'.\n")

def test_baltica_cli_absent_config():
    runner = CliRunner()
    result = runner.invoke(cli, 'rmats NOCONFIG')
    assert result.exit_code == 2
    assert result.output.endswith("Error: Invalid value for 'CONFIG_FILE': Path 'NOCONFIG' does not exist.\n")

def test_avaiable_workflows_none(fs):
    baltica_path = Path('.')
    fs.create_file('a.tmp')
    assert avaiable_workflows(baltica_path) == []

def test_avaiable_workflows_single(fs):
    baltica_path = Path('.')
    fs.create_file('a.smk')
    assert avaiable_workflows(baltica_path) == ['a']

def test_avaiable_workflows_multiple(fs):
    baltica_path = Path('.')
    fs.create_file('a.smk')
    fs.create_file('b.smk')
    fs.create_file('c.smk')
    assert avaiable_workflows(baltica_path) == ['a', 'b', 'c']


# def test
# mock baltica installation to be withint test
# only workflow avaiable  should be  test.smk
# snakemake -s .tests/test_1.smk -j1 --use-singularity
# should work with and without and test singularity
# How to test something only avaiable withint singularity? godlovedc/lolcow !!!
# 
# subprocess.checkoutput
# @patch("baltica.command.subprocess.run")
# def test_singularity_unavaiable(mock_run):
#     mock_stdout = MagicMock()
#     mock_stdout.configure_mock(
#         **{
#             "stdout.decode.return_value": '{"A": 3}'
#         }
#     )
#     mock_run.return_value = mock_stdout
#     runner = CliRunner()
#     result = runner.invoke(cli, '--help')
# @patch("baltica.command.subprocess.run")
# def test_singularity_avaiable():
#     pass
#     result = runner.invoke(cli, 'rmats')
