# -*- coding: utf-8 -*-
#
#  This file is part of Bioconvert software
#
#  Copyright (c) 2017 - Bioconvert Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
#
##############################################################################
""" description """
import json
import glob
from os.path import join, basename, exists
from easydev import md5
import os


def download_singularity_image(outfile, container_path, md5value=None, force=False):
    assert outfile.endswith(".simg"), "output filename must be .simg"

    # download singularity
    from bioconvert import configuration as config
    # note that in singularity v2.4, whatever extension you put, it is
    # replaced by simg
    singfile = "{}/{}".format(config.user_config_dir, outfile)

    if exists(singfile) and md5value and md5(singfile) == md5value and force is False:
        print("Found singularity (graphviz) image")
    elif exists(singfile) and md5value is None and force is False:
        print("Found singularity (graphviz) image but md5 not checked")
    else:
        print("Downloading singularity. Please wait")
        cmd = "singularity pull --name {}  {}"
        if force is True:
            cmd += " -F "
        cmd = cmd.format(singfile, container_path)
        print(cmd)
        try:
            shell(cmd)
        except:
            import os
            os.system(cmd)
    return singfile














