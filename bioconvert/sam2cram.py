"""Convert :term:`SAM` file to :term:`CRAM` file"""
import os

from bioconvert import ConvBase

import colorlog
logger = colorlog.getLogger(__name__)

__all__ = ["SAM2CRAM"]

class SAM2CRAM(ConvBase):
    """Convert :term:`SAM` file to :term:`CRAM` file

    The conversion requires the reference corresponding to the input file
    It can be provided as an argument in the constructor. Otherwise, 
    a local file with same name as the input file but an .fa extension is looked
    for. Otherwise, we ask for the user to provide the input file. This is 
    useful for the standalone application.

    """
    input_ext = [".sam"]
    output_ext = [".cram"]

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

        self._default_method = "samtools"

        self.reference = reference
        if self.reference is None:
            logger.debug("No reference provided. Infering from input file")
            # try to find the local file replacing .sam by .fa
            reference = infile.replace(".sam", ".fa")
            if os.path.exists(reference):
                logger.debug("Reference found from inference ({})".format(reference))
            else:
                logger.debug("No reference found.")
                msg = "Please enter the reference corresponding "
                msg += "to the input SAM file:"
                reference = input(msg)
                if os.path.exists(reference) is False:
                    raise IOError("Reference required")
                else:
                    logger.debug("Reference exist ({}).".format(reference))

            self.reference = reference

    def _method_samtools(self, *args, **kwargs):
        # -S means ignored (input format is auto-detected)
        # -b means output to BAM format
        # -h means include header in SAM output
        cmd = "samtools view -@ {} -SCh {} -T {} > {}".format(self.threads, 
            self.infile, self.reference, self.outfile)
        try:
            self.execute(cmd)
        except:
            logger.debug("FIXME. The ouput message from samtools is on stderr...")
