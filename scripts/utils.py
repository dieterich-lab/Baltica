#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 11:13 2019-03-14 2019 

"""
from pathlib import Path


def basename(path, suffix=None):
    if suffix:
        return str(Path(path).with_suffix(suffix).name)
    return str(Path(path).name)
