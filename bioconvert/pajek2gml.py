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
"""Convert :term:`PAJEK` to :term:`GML` format"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires

__all__ = ["PAJEK2GML"]


class PAJEK2GML(ConvBase):
    """Convert :term:`PAJEK` file into :term:`GML` file

    Conversion based on the :class:`networkx` library.

    .. plot::

        from bioconvert.pajek2gml import PAJEK2GML
        img = PAJEK2GML.example_bench(N=5)
        img.savefig("bench_pajek2gml.png", dpi=200, bbox_inches="tight")
    """

    #: Default value
    _default_method = "networkx"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input :term:`PAJEK` file
        :param str outfile: output :term:`GML` file
        """
        super(PAJEK2GML, self).__init__(infile, outfile, *args, **kargs)

    @requires(python_library="networkx")
    def _method_networkx(self, *args, **kwargs):
        """Convert :term:`PAJEK` to :term:`GML` using :class:`networkx`"""
        import networkx as nx

        G = nx.read_pajek(self.infile)
        nx.write_gml(G, self.outfile)
