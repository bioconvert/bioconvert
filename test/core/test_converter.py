from bioconvert import Bioconvert
from easydev import TempFile

import pytest

from bioconvert.bam2fasta import BAM2FASTA
from bioconvert.gz2bz2 import GZ2BZ2

from .. import test_dir

@pytest.mark.skipif(len(BAM2FASTA.available_methods) == 0, reason="missing dependencies")
def test_bioconvert():
    infile = f"{test_dir}/data/bam/test_measles.sorted.bam"
    with TempFile(suffix=".fasta") as fout:
        c = Bioconvert(infile, fout.name, force=True)
        c()


def test_bioconvert_force_false():
    infile = f"{test_dir}/data/fastq/ERR3295124.fastq"
    with TempFile(suffix=".fasta") as fout:
        c = Bioconvert(infile, fout.name, force=True)
        c()
        with pytest.raises(ValueError):
            Bioconvert(infile, fout.name, force=False)


def test_bioconvert_dsrc_only_for_fastq():
    with TempFile(suffix=".fastq.dsrc") as fout, TempFile(suffix=".fasta") as fin:
        c = Bioconvert(fin.name, fout.name, force=True)
        assert c is not None
    with TempFile(suffix=".fasta.dsrc") as fout, TempFile(suffix=".fasta") as fin:
        with pytest.raises(IOError):
            Bioconvert(fin.name, fout.name, force=True)


@pytest.mark.skipif(len(GZ2BZ2.available_methods) == 0, reason="missing dependencies")
def test_bioconvert_decompression_compression_mode():
    infile = f"{test_dir}/data/gz/measles_R1.fastq.gz"
    with TempFile(suffix=".fastq.bz2") as fout:
        c = Bioconvert(infile, fout.name, force=True)


def test_indirect_conversion_impossible():
    with TempFile(suffix=".json") as fout, TempFile(suffix=".paf") as fin:
        with pytest.raises(Exception) as err:
            Bioconvert(fin.name, fout.name, force=True,)
            assert (str(err.value)
                    == "Requested input format ('JSON') to output format ('PAF') is not available in bioconvert")


def test_indirect_conversion():
    infile = f"{test_dir}/data/fastq/ERR3295124.fastq"
    with TempFile(suffix=".clustal") as fout:
        c = Bioconvert(infile, fout.name, force=True)
        c()