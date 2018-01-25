import pytest
from bioconvert.core.registry import Registry


def test_registry():
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
