from pathlib import Path
import shlex
import subprocess
import sys

import pytest


@pytest.mark.parametrize(
    ("strandness", "expected_flag"),
    [
        ("fr-firststrand", "--rf"),
        ("fr-secondstrand", "--fr"),
        ("unstranded", ""),
        (None, ""),
    ],
)
def test_stringtie_dry_run_uses_configured_strand_flag(
    tmp_path, strandness, expected_flag
):
    """StringTie receives the matching flag for the configured library type."""
    workflow = Path(__file__).parents[1] / "baltica" / "stringtie.smk"
    (tmp_path / "stringtie" / "merged_bam").mkdir(parents=True)
    (tmp_path / "stringtie" / "merged_bam" / "condition.bam").touch()

    config = [
        "samples:",
        "  condition_rep: sample.bam",
        "sample_path: samples/",
        "ref: reference.gtf",
    ]
    if strandness is not None:
        config.append(f"strandness: {strandness}")
    config_path = tmp_path / "config.yaml"
    config_path.write_text("\n".join(config) + "\n")

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--snakefile",
            str(workflow),
            "--directory",
            str(tmp_path),
            "--configfile",
            str(config_path),
            "--cores",
            "1",
            "--dry-run",
            "--printshellcmds",
            "stringtie/stringtie/condition.gtf",
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    commands = [
        line
        for line in (result.stdout + result.stderr).splitlines()
        if "stringtie stringtie/merged_bam/condition.bam" in line
    ]
    assert len(commands) == 1
    strand_flags = {"--rf", "--fr"} & set(shlex.split(commands[0]))
    assert strand_flags == ({expected_flag} if expected_flag else set())
