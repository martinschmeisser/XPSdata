# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 20:55:56 2018

@author: jdo
"""

from XPS.SpecsXMLFile import SpecsXMLFile
from XPS.SpecsSLEFile import SpecsSLEFile


def get_region(filename, group_no, region_no):
    """ returns a region class object from a Specs XPS data file

    load from .xml (SpecsLab2) or .sle (SpecsLab Prodigy) files"""

    if filename[-3:].lower() == 'xml':
        data = SpecsXMLFile(filename, encoding='ansi')
    elif filename[-3:].lower() == 'sle':
        data = SpecsSLEFile(filename)

    return data[group_no][region_no]
