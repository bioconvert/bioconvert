
from bioconvert.core.registry import Registry



def test_registry():
    rr = Registry()
    assert rr.conversion_exists("bam", "bed") is True
    assert rr.conversion_exists("bam", "dummy") is False
