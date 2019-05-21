# -*- coding: utf-8 -*-
###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Authors: see CONTRIBUTORS.rst                                           #
# Copyright © 2018  Institut Pasteur, Paris and CNRS.                     #
# See the COPYRIGHT file for details                                      #
#                                                                         #
# bioconvert is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# bioconvert is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program (COPYING file).                                 #
# If not, see <http://www.gnu.org/licenses/>.                             #
###########################################################################
"""Download singularity image"""
from os.path import exists
from easydev import md5


__all__ = ['download_singularity_image']


def download_singularity_image(outfile, container_path, md5value=None,
                               force=False):

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
            from bioconvert.core.shell import shell
            shell(cmd)
        except:
            import os
            os.system(cmd)
    return singfile














