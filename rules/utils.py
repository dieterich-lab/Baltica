import re

from snakemake.io import expand




def natural_sort_key(s, _nsre=re.compile("([0-9]+)")):
    return [int(text) if text.isdigit() else text.lower()
            for text in _nsre.split(s)]


def comparison(wc, index, mapping):
    condition = wc.contrast.split("-vs-")[index]
    return expand("majiq/{name}.majiq", name=mapping[condition])
