import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import requires
from bioconvert.core.utils import generate_outfile_name

_log = colorlog.getLogger(__name__)


class BPLINK2PLINK(ConvBase):
    """Converts a genotype dataset bed+bim+fam in :term:`BPLINK` format to
    ped+map :term:`PLINK` format

    Conversion is based on plink executable

    """
    _default_method = 'plink'

    def __init__(self, infile, outfile=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`BPLINK` file.
        :param str outfile: (optional) output :term:`PLINK` file
        """
        if not outfile:
            outfile = generate_outfile_name(infile, '')
        super().__init__(infile, outfile)

    @requires("plink")
    def _method_plink(self, *args, **kwargs):
        """
        Convert plink file in text using plink executable.
        """
        cmd = 'plink --bfile {infile} --recode --out {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
