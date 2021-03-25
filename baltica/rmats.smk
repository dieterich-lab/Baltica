# -*- coding: utf-8
"""
Created on 17:07 27/07/2018
Snakemake workflow for rMATS
TODO citation
TODO notes 
"""
__author__ = "Thiago Britto Borges"
__copyright__ = "Copyright 2020, Dieterichlab"
__email__ = "Thiago.BrittoBorges@uni-heidelberg.de"
__license__ = "MIT"

from collections import defaultdict
import tempfile


workdir: config.get("path", ".")
contrasts = config['contrasts']
keys = config['samples'].keys()
keys = [tuple(x.split('_')) for x in keys]

shell.prefix('source ~/miniconda3/etc/profile.d/conda.sh; conda init; conda activate rmats;')
temp_dir = tempfile.TemporaryDirectory()
strand = {
  'reverse': 'fr-secondstrand',
  'forward': 'fr-firststrand'}

d = defaultdict(list)
for x in keys: d[x[0]].append(x[1])

include: "symlink.smk"


rule all:
    input:
    	expand("rmats/{group}.txt", group=d.keys()),
        expand("rmats/{contrast}/", contrast=contrasts)


rule create_rmats_input:
  input: lambda wc: expand("mappings/{{group}}_{rep}.bam", rep=d.get(wc.group))
  output: "rmats/{group}.txt"
  run:
  	with open(str(output), "w") as fou:
  		fou.write(','.join(input))
                #fou.write(config['read_len'])


rule run_rmats:
  input: "rmats/{alt}.txt", "rmats/{ref}.txt",
  output: directory("rmats/{alt}-vs-{ref}/")
  shadow: "shallow"
  threads: 10
  params:
      gtf = config['ref'],
      # is_paired = '-t single' if config['is_single_end'] else '',
      lib = "--libType " + strand.get(config['strandness'], 'fr-unstranded'),
      read_len = config['read_len'],
      allow_clipping = "--allow-clipping",
      tmp = os.path.join(temp_dir.name, '{alt}_vs_{ref}/')
  shell: 
      "python ~/rmats-turbo/rmats.py "
      "--b1 {input[0]} "
      "--b2 {input[1]} "
      "--gtf {params.gtf} "
      "--variable-read-length " 
      "--readLength {params.read_len} "
      "--nthread {threads} " 
      "--novelSS "
      "{params.lib} " 
      "{params.allow_clipping} "
      "--od {output} " 
      "--tmp {params.tmp}" 
