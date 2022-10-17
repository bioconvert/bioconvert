###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Copyright © 2018-2022  Institut Pasteur, Paris and CNRS.                #
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
#                                                                         #
# Repository: https://github.com/bioconvert/bioconvert                    #
# Documentation: http://bioconvert.readthedocs.io                         #
###########################################################################
""" Convert a compressed :term:`FASTQ` from :term:`DSRC` to :term:`FASTQ` format"""
import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

logger = colorlog.getLogger(__name__)


__all__ = ["DSRC2GZ"]


class DSRC2GZ(ConvBase):
    """Convert a compressed :term:`FASTQ` from :term:`DSRC` to :term:`GZ` format

    Methods available are based on dsrc [DSRC]_ and pigz [PIGZ]_.

    """

    #: Default value
    _default_method = "dsrcpigz"
    _threading = True

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input DSRC filename
        :param str outfile: output GZ filename

        """
        super(DSRC2GZ, self).__init__(infile, outfile, *args, **kargs)

    @requires("dsrc")
    def _method_dsrcpigz(self, *args, **kwargs):
        """Do the conversion dsrc -> :term:`GZ`.
        Method that uses pigz and dsrc.

        `pigz documentation <https://linux.die.net/man/1/pigz>`_
        `dsrc documentation <https://github.com/refresh-bio/DSRC>`_

        option threadig does not work with the dsrc version from conda so we
        do not add the -t threads option
        """
        cmd = "dsrc d -s  {input} | pigz -c -p {threads} > {output}"
        self.execute(cmd.format(threads=self.threads, input=self.infile, output=self.outfile))
