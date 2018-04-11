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

"""Convert :term:`SAM` file to :term:`CRAM` file"""
import os

from bioconvert import ConvBase

import colorlog

from bioconvert.core.decorators import requires

logger = colorlog.getLogger(__name__)

__all__ = ["SAM2CRAM"]


class SAM2CRAM(ConvBase):
    """Convert :term:`SAM` file to :term:`CRAM` file

    The conversion requires the reference corresponding to the input file
    It can be provided as an argument in the constructor. Otherwise,
    a local file with same name as the input file but an .fa extension is
    looked for. Otherwise, we ask for the user to provide the input file.
    This is useful for the standalone application.

    """
    _default_method = "samtools"

    def __init__(self, infile, outfile, reference=None, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input SAM file
        :param str outfile: output filename
        :param str reference: reference file in :term:`FASTA` format

        command used::

            samtools view -SCh

        .. note:: the API related to the third argument may change in the future.
        """
        super(SAM2CRAM, self).__init__(infile, outfile, *args, **kargs)


        self.reference = reference
        if self.reference is None:
            logger.warning(
                "No reference provided. Looking for same file "
                "with .fa extension in same directory")
            # try to find the local file replacing .sam by .fa
            reference = infile.replace(".sam", ".fa")
            if os.path.exists(reference):
                logger.info(
                    "Reference found from inference ({})".format(reference))
            else:
                logger.warning("No reference found.")
                msg = "Please enter the reference corresponding "
                msg += "to the input SAM file:"
                reference = input(msg)
                if os.path.exists(reference) is False:
                    raise IOError("Reference required")
                else:
                    logger.debug("Reference exists ({}).".format(reference))

            self.reference = reference

    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        # -C means output is CRAM
        reference = kwargs.get("reference", self.reference)

        cmd = "samtools view -@ {} -C {} -T {} > {}".format(
            self.threads, self.infile, reference, self.outfile)
        try:
            self.execute(cmd)
        except:
            logger.debug("FIXME. The ouput message from samtools is on stderr...")
