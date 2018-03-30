import pytest

from bioconvert.bam2bed import BAM2BED
from bioconvert.core.registry import Registry
from bioconvert.sra2fastq import SRA2FASTQ


@pytest.mark.skipif(len(BAM2BED.available_methods) == 0
                    or len(SRA2FASTQ.available_methods) == 0,
                    reason="missing dependencies")
def test_registry_with_dependencies():
    rr = Registry()
    # Isnt'it supposed to be upper case?
    # assert rr.conversion_exists("bam", "bed") is True
    # assert rr.conversion_exists("bam", "dummy") is False
    assert rr.conversion_exists("BAM", "BED")
    assert not rr.conversion_exists("BAM", "DUMMY")
    assert rr.conversion_exists("SRA", "BAM", allow_indirect=True)

    assert ('BAM', 'BED') in rr
    assert ('NIMPORT', 'NAOIK') not in rr

    converter = rr[('BAM', 'BED')]
    with pytest.raises(KeyError) as err:
        rr[('BAM', 'BED')] = converter
    assert (str(err.value)
            == "'an other converter already exists for BAM -> BED'")


def test_registry_with_less_dependencies():
    rr = Registry()
    # Isnt'it supposed to be upper case?
    # assert rr.conversion_exists("bam", "bed") is True
    # assert rr.conversion_exists("bam", "dummy") is False
    assert rr.conversion_exists("FASTQ", "FASTA")
    assert not rr.conversion_exists("BAM", "DUMMY")
    assert rr.conversion_exists("CLUSTAL", "FASTQ", allow_indirect=True)

    assert ('FASTQ', 'FASTA') in rr
    assert ('NIMPORT', 'NAOIK') not in rr

    converter = rr[('FASTQ', 'FASTA')]
    with pytest.raises(KeyError) as err:
        rr[('FASTQ', 'FASTA')] = converter
    assert (str(err.value)
            == "'an other converter already exists for FASTQ -> FASTA'")
