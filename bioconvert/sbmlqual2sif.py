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
"""Convert :term:`SBMLQUAL` format to :term:`SIF` format"""
import colorlog

from bioconvert.core.base import ConvBase
from bioconvert.core.decorators import requires

logger = colorlog.getLogger(__name__)

__all__ = ["SBMLQUAL2SIF"]


class SBMLQUAL2SIF(ConvBase):
    """Convert :term:`SBMLQUAL` file into :term:`SIF` file

    :term:`SBML` qual is an SBML Level 3 package for qualitative models
    (Boolean networks, logical regulatory models). :term:`SIF` (Simple
    Interaction Format) is used by Cytoscape to represent networks.

    The conversion extracts pairwise interactions from SBML qual transitions
    and writes them as ``source<tab>interaction<tab>target`` triples.  The
    interaction sign is derived from the ``sign`` attribute of each
    :class:`libsbml.Input`:

    * ``positive``  → ``activates``
    * ``negative``  → ``inhibits``
    * ``dual``      → ``regulates``
    * unset / other → ``regulates``

    Species that appear in no transition are also written as isolated nodes.

    Available methods: python_libsbml

    .. seealso:: `SBML qual specification
        <http://www.sbml.org/Documents/Specifications/SBML_Level_3/Packages/qual>`_
    """

    #: Default value
    _default_method = "python_libsbml"

    def __init__(self, infile, outfile):
        """.. rubric:: Constructor

        :param str infile: input :term:`SBMLQUAL` file.
        :param str outfile: output :term:`SIF` file.
        """
        super(SBMLQUAL2SIF, self).__init__(infile, outfile)

    @requires(python_library="python-libsbml")
    def _method_python_libsbml(self, *args, **kwargs):
        """Convert :term:`SBMLQUAL` to :term:`SIF` using *python-libsbml*.

        `python-libsbml documentation
        <https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/>`_
        """
        import libsbml

        _SIGN_MAP = {
            libsbml.INPUT_SIGN_POSITIVE: "activates",
            libsbml.INPUT_SIGN_NEGATIVE: "inhibits",
            libsbml.INPUT_SIGN_DUAL: "regulates",
        }

        doc = libsbml.readSBMLFromFile(self.infile)
        errors = doc.getNumErrors()
        for i in range(errors):
            err = doc.getError(i)
            if err.getSeverity() >= libsbml.LIBSBML_SEV_ERROR:
                raise ValueError(
                    "SBML read error in {}: {}".format(self.infile, err.getMessage())
                )

        model = doc.getModel()
        if model is None:
            raise ValueError("No model found in SBML file: {}".format(self.infile))

        mplugin = model.getPlugin("qual")
        if mplugin is None:
            raise ValueError(
                "No 'qual' package found in SBML model. "
                "Make sure the file is a valid SBML qual document."
            )

        # Collect all species ids for isolated-node tracking
        all_species = set()
        for i in range(mplugin.getNumQualitativeSpecies()):
            qs = mplugin.getQualitativeSpecies(i)
            all_species.add(qs.getId())

        written_species = set()
        interactions = []

        for i in range(mplugin.getNumTransitions()):
            transition = mplugin.getTransition(i)
            num_inputs = transition.getNumInputs()
            num_outputs = transition.getNumOutputs()

            for j in range(num_inputs):
                inp = transition.getInput(j)
                source = inp.getQualitativeSpecies()
                sign = inp.getSign()
                interaction = _SIGN_MAP.get(sign, "regulates")

                for k in range(num_outputs):
                    out = transition.getOutput(k)
                    target = out.getQualitativeSpecies()
                    interactions.append((source, interaction, target))
                    written_species.add(source)
                    written_species.add(target)

        with open(self.outfile, "w") as fh:
            for source, interaction, target in interactions:
                fh.write("{}\t{}\t{}\n".format(source, interaction, target))
            # Write isolated species (not part of any transition)
            for species in sorted(all_species - written_species):
                fh.write("{}\n".format(species))
