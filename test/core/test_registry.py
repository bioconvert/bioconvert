import pytest
from bioconvert.core.registry import Registry



def test_registry():
    rr = Registry()
    assert rr.conversion_exists("bam", "bed") is True
    assert rr.conversion_exists("bam", "dummy") is False

    assert ('BAM', 'BED') in rr
    assert ('NIMPORT', 'NAOIK') not in rr

    converter = rr[('BAM', 'BED')]
    with pytest.raises(KeyError) as err:
        rr[('BAM', 'BED')] = converter
    assert "'an other converter already exist for BAM -> BED'" == str(err.value)
