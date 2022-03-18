import pytest

from bioconvert.bam2fasta import BAM2FASTA
from easydev import TempFile, md5

from . import test_dir

#@pytest.mark.skipif(BAM2FASTA._method_samtools.is_disabled, reason="missing dependencies")
@pytest.mark.parametrize("method", BAM2FASTA.available_methods)
def test_methods(method):
    infile = f"{test_dir}/data/bam/test_measles.sorted.bam"

    with TempFile(suffix=".fa") as tempfile:
        convert = BAM2FASTA(infile, tempfile.name)
        convert(method=method)
        # samtools 1.6 / hstlib 1.6 gives different results on travis and
        # locally
        assert md5(tempfile.name.replace(".", "_1.")) in [
            "9242d127969a089ddeedbc2002c62686"]
        assert md5(tempfile.name.replace(".", "_2.")) in [
            "b753ad368c9614130884acb29861bd23"]

    for ext in ['gz', 'bz2']:# Test compression 
        with TempFile(suffix=".fasta.{}".format(ext)) as tempfile:
            convert = BAM2FASTA(infile, tempfile.name)
            convert(method=method)
            # no check, just running 

    infile = f"{test_dir}/data/bam/test_measles_unpaired.sorted.bam"

    with TempFile(suffix=".fa") as tempfile:
        convert = BAM2FASTA(infile, tempfile.name)
        convert(method=method)
        # samtools 1.6 / hstlib 1.6 gives different results on travis and
        # locally

    for ext in ['gz', 'bz2']:# Test compression 
        with TempFile(suffix=".fasta.{}".format(ext)) as tempfile:
            convert = BAM2FASTA(infile, tempfile.name)
            convert(method=method)
            # no check, just running 