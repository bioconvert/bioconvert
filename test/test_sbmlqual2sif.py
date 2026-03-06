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
import pytest

from bioconvert import TempFile, md5
from bioconvert.sbmlqual2sif import SBMLQUAL2SIF

from . import test_dir


@pytest.mark.parametrize("method", SBMLQUAL2SIF.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/sbmlqual/test_model.sbmlqual"
    expected_outfile = f"{test_dir}/data/sif/test_model.sif"
    with TempFile(suffix=".sif") as tempfile:
        converter = SBMLQUAL2SIF(infile, tempfile.name)
        converter(method=method)
        assert md5(tempfile.name) == md5(expected_outfile)
