import pytest

from bioconvert.bam2cov import BAM2COV
from bioconvert.core.registry import Registry
from bioconvert.sra2fastq import SRA2FASTQ


@pytest.mark.skipif(len(BAM2COV.available_methods) == 0
                    or len(SRA2FASTQ.available_methods) == 0,
                    reason="missing dependencies")
def test_registry_with_dependencies():
    rr = Registry()
    assert rr.conversion_exists(("BAM",), ("COV",))
    assert not rr.conversion_exists(("BAM",), ("DUMMY",))
    assert rr.conversion_exists(("SRA",), ("BAM",), allow_indirect=True)

    assert (('BAM',), ('COV',)) in rr
    assert (('NIMPORT',), ('NAOIK',)) not in rr

    converter = rr[('BAM',), ('COV',)]
    with pytest.raises(KeyError) as err:
        rr[('BAM',), ('COV',)] = converter
    assert (str(err.value)
            == "'an other converter already exists for BAM -> COV'")


def test_registry_with_less_dependencies():
    rr = Registry()
    assert rr.conversion_exists(("FASTQ",), ("FASTA",))
    assert not rr.conversion_exists(("BAM",), ("DUMMY",))
    assert rr.conversion_exists(("CLUSTAL",), ("FASTQ",), allow_indirect=True)

    # one to one
    assert (('FASTQ',), ('FASTA',)) in rr
    # unexisting entry
    assert (('NIMPORT',), ('NAOIK',)) not in rr
    # one to many
    assert (('FASTQ',), ('FASTA', 'QUAL')) in rr
    # without tuples
    assert ['FASTQ', 'FASTA'] in rr
    assert [('FASTQ',), 'FASTA'] in rr
    assert [['FASTQ'], 'FASTA'] in rr
    assert [['FASTQ'], ['FASTA']] in rr
    assert ['FASTQ', ('FASTA',)] in rr
    # argument must be a list or tuple
    try:
        'FASTA' in rr
        assert False
    except ValueError:
        assert True
    # must have 2 items only
    try:
        ['a', 'b', 'c'] in rr
        assert False
    except ValueError:
        assert True

    # iterator through the conversions
    for this in rr:
        assert this

    converter = rr[('FASTQ',), ('FASTA',)]
    with pytest.raises(KeyError) as err:
        rr[('FASTQ',), ('FASTA',)] = converter
    assert (str(err.value)
            == "'an other converter already exists for FASTQ -> FASTA'")


def test_rgistry():
    rr = Registry()
    rr.info()
    print(rr)





