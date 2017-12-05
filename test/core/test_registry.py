import pytest
from bioconvert.core.registry import Registry



def test_registry():
    rr = Registry()
    assert rr.conversion_exists("bam", "bed") is True
    assert rr.conversion_exists("bam", "dummy") is False

    assert ('.bam', '.bed') in rr
    assert ('.nimport', '.naoik') not in rr

    converter = rr[('.bam', '.bed')]
    with pytest.raises(KeyError) as err:
        rr[('.bam', '.bed')] = converter
    assert "'an other converter already exist for .bam -> .bed'" == str(err.value)
