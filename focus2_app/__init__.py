# -*- coding: utf-8 -*-

import pkg_resources

try:
    version = pkg_resources.require("metagenomics-focus2")[0].version
except:
    raise ValueError('Cannot find version number')
