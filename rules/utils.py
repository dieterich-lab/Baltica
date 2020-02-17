import re

from snakemake.io import expand


def extract_samples_replicates(samples, _pattern=re.compile('^(.+)_(.+)$')):
    """
    Extract pairs of condition and replicate name from sample files

    :param str _pattern: pattern to . Default uses {condition}_{replicate} template
    :param list samples:
    :return:
    :rtype: list
    """
    return list(zip(*[re.match(_pattern, x).groups() for x in samples]))


def rename_bam_to_fastq(x):
    # reads star log for input files
    f = x.replace("Aligned.noS.bam", "Log.out")
    # f.replace(f.name, 'Log.out'
    with open(f) as fin:
        for line in fin:
            if line.startswith("##### Command Line:"):
                line = next(fin)
                break

    for args in line.split("--"):
        if args.startswith("readFilesIn"):
            break
    return [x.split("mapping")[0] + arg for arg in args.split()[1:3]]


def natural_sort_key(s, _nsre=re.compile("([0-9]+)")):
    return [int(text) if text.isdigit() else text.lower()
            for text in _nsre.split(s)]


def comparison(wc, index, mapping):
    condition = wc.contrast.split("-vs-")[index]
    return expand("majiq/{name}.majiq", name=mapping[condition])
