"""Convert :term:`SAM` file to :term:`CRAM` file"""
import os
from bioconvert import ConvBase

import colorlog
_log = colorlog.getLogger(__name__)


class CRAM2SAM(ConvBase):
    """Convert :term:`CRAM` file to :term:`SAM` file

    The conversion requires the reference corresponding to the input file
    It can be provided as an argument in the constructor. Otherwise, 
    a local file with same name as the input file but an .fa extension is looked
    for. Otherwise, we ask for the user to provide the input file. This is 
    useful for the standalone application.

    """
    input_ext = [".cram"]
    output_ext = ".sam"

    def __init__(self, infile, outfile, reference=None, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input SAM file
        :param str outfile: output filename
        :param str reference: reference file in :term:`FASTA` format

        command used::

            samtools view -@ <thread> -Sh -T <reference> in.cram > out.sam

        .. note:: the API related to the third argument may change in the future.
        """
        super(self).__init__(infile, outfile, *args, **kargs)

        self._default_method = "samtools"
        self.reference = reference

        if self.reference is None:
            _log.debug("No reference provided. Infering from input file")
            # try to find the local file replacing .sam by .fa
            reference = infile.replace(".cram", ".fa")
            if os.path.exists(reference):
                _log.debug("Reference found from inference ({})".format(reference))
            else:
                _log.debug("No reference found.")
                msg = "Please enter the reference corresponding "
                msg += "to the input SAM file:"
                reference = input(msg)
                if not os.path.exists(reference):
                    raise IOError("Reference required")
                else:
                    _log.debug("Reference exist ({}).".format(reference))

            self.reference = reference


    def _method_samtools(self, threads=None):
        # -S means ignored (input format is auto-detected)
        # -h means include header in SAM output
        if threads is None:
            threads = self.max_threads
        cmd = "samtools view -@ {thread} -Sh -T {ref} {infile} > {outfile}".format(
            thread=threads,
            ref=self.reference,
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)





