sample_path = config["sample_path"]
name = config["samples"].keys()
sample = config["samples"].values()

from os import symlink, path
from shutil import copy2

def make_symlink(origin, dest):
    relative_path = path.relpath(
        path.dirname(origin), 
        path.dirname(dest)
    )
    new_origin = path.join(relative_path, path.basename(dest))

    symlink(new_origin, dest)

rule symlink:
    input:
        bam=expand("{samples_dir}{sample}", sample=sample, samples_dir=sample_path),
        bai=expand("{samples_dir}{sample}.bai", sample=sample, samples_dir=sample_path)
    output:
        bam=expand('mappings/{name}.bam', name=name),
        bai=expand('mappings/{name}.bam.bai', name=name)
    run:
        for bam_in, bai_in, bam_out, bai_out in zip(
                input.bam, input.bai, output.bam, output.bai):
            copy2(bam_in, bam_out)
            copy2(bai_in, bai_out)
