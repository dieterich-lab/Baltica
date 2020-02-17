import os

sample_path = config["sample_path"]

rule symlink:
    input: bam=expand("{samples_dir}{name}", name=name, samples_dir=sample_path),
         bai=expand("{samples_dir}{name}.bai", name=name, samples_dir=sample_path)
    output: bam=expand('mappings/{sample}.bam', sample=sample),
          bai=expand('mappings/{sample}.bam.bai', sample=sample)
    run:
        for bam_in, bai_in, bam_out, bai_out in zip(
                input.bam, input.bai, output.bam, output.bai):
            os.symlink(bam_in, bam_out)
            os.symlink(bai_in, bai_out)

rule create_log:
    output: directory('logs')
    run:
        try:
            os.mkdir('logs/')
        except os.OSError.FileExistsError:
            pass
