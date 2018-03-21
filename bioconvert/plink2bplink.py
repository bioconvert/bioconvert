import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import requires
from bioconvert.core.utils import generate_outfile_name

_log = colorlog.getLogger(__name__)


class PLINK2BPLINK(ConvBase):
    """Converts a genotype dataset ped+map in :term:`PLINK` format to
    bed+bim+fam :term:`BPLINK` format

    Conversion is based on plink executable

    """
    _default_method = 'plink'

    def __init__(self, infile, outfile=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`PLINK` file.
        :param str outfile: (optional) output :term:`BPLINK` file
        """
        if not outfile:
            outfile = generate_outfile_name(infile, '')
        super().__init__(infile, outfile)

    @requires("plink")
    def _method_plink(self, *args, **kwargs):
        """
        Convert plink file in text using plink executable.
        """
        cmd = 'plink --file {infile} --make-bed --out {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
