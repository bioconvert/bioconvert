import re
import pytest
import os

from bioconvert.core.decorators import requires, get_known_dependencies_with_availability

skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ, reason="On travis")


### AUTOMATICALLY GENERATED TESTS (START)

def test_require_awk():
    assert requires(external_binary="awk")(object()).is_disabled is False


def test_require_bamCoverage():
    assert requires(external_binary="bamCoverage")(object()).is_disabled is False


def test_require_bamtools():
    assert requires(external_binary="bamtools")(object()).is_disabled is False


def test_require_bcftools():
    assert requires(external_binary="bcftools")(object()).is_disabled is False


def test_require_bedGraphToBigWig():
    assert requires(external_binary="bedGraphToBigWig")(object()).is_disabled is False


def test_require_bedtools():
    assert requires(external_binary="bedtools")(object()).is_disabled is False


def test_require_bigWigToBedGraph():
    assert requires(external_binary="bigWigToBedGraph")(object()).is_disabled is False


def test_require_bioawk():
    assert requires(external_binary="bioawk")(object()).is_disabled is False


def test_require_bunzip2():
    assert requires(external_binary="bunzip2")(object()).is_disabled is False


def test_require_conda():
    assert requires(external_binary="conda")(object()).is_disabled is False


def test_require_dsrc():
    assert requires(external_binary="dsrc")(object()).is_disabled is False


def test_require_faToTwoBit():
    assert requires(external_binary="faToTwoBit")(object()).is_disabled is False


def test_require_fastq_dump():
    assert requires(external_binary="fastq-dump")(object()).is_disabled is False


def test_require_ls():
    assert requires(external_binary="ls")(object()).is_disabled is False


def test_require_mawk():
    assert requires(external_binary="mawk")(object()).is_disabled is False


def test_require_mv():
    assert requires(external_binary="mv")(object()).is_disabled is False


def test_require_perl():
    assert requires(external_binary="perl")(object()).is_disabled is False


def test_require_pigz():
    assert requires(external_binary="pigz")(object()).is_disabled is False


@skiptravis
def test_require_plink():
    assert requires(external_binary="plink")(object()).is_disabled is False


def test_require_rm():
    assert requires(external_binary="rm")(object()).is_disabled is False


def test_require_sambamba():
    assert requires(external_binary="sambamba")(object()).is_disabled is False


def test_require_samtools():
    assert requires(external_binary="samtools")(object()).is_disabled is False


def test_require_sed():
    assert requires(external_binary="sed")(object()).is_disabled is False


def test_require_seqtk():
    assert requires(external_binary="seqtk")(object()).is_disabled is False


def test_require_squizz():
    assert requires(external_binary="squizz")(object()).is_disabled is False


def test_require_top():
    assert requires(external_binary="top")(object()).is_disabled is False


def test_require_twoBitToFa():
    assert requires(external_binary="twoBitToFa")(object()).is_disabled is False


def test_require_Bio():
    assert requires(python_library="Bio")(object()).is_disabled is False


def test_require_csv():
    assert requires(python_library="csv")(object()).is_disabled is False


def test_require_gatb():
    assert requires(python_library="gatb")(object()).is_disabled is False


def test_require_mappy():
    assert requires(python_library="mappy")(object()).is_disabled is False


def test_require_mmap():
    assert requires(python_library="mmap")(object()).is_disabled is False


def test_require_os():
    assert requires(python_library="os")(object()).is_disabled is False


def test_require_pandas():
    assert requires(python_library="pandas")(object()).is_disabled is False


def test_require_pysam():
    assert requires(python_library="pysam")(object()).is_disabled is False


def test_require_sys():
    assert requires(python_library="sys")(object()).is_disabled is False


### AUTOMATICALLY GENERATED TESTS (END)


def test_require_all_and_print_test():
    known_missing_dependencies = ["tagada%i" % i for i in range(10)]
    known_missing_dependencies.append("plink")

    assert type(get_known_dependencies_with_availability(as_dict=True)) == dict

    print("""
### AUTOMATICALLY GENERATED TESTS (START)
""")
    for d, s, t in get_known_dependencies_with_availability():
        if d[0:6] != "tagada":
            if d in known_missing_dependencies:
                print("@skiptravis")
            print("""def test_require_{}():
    assert requires({}="{}")(object()).is_disabled is False

""".format(re.sub('[^a-zA-Z0-9]', "_", d), "external_binary" if t == "binary" else "python_library", d))
    print("""### AUTOMATICALLY GENERATED TESTS (END)
""")
    for d, s, t in get_known_dependencies_with_availability():
        assert d in known_missing_dependencies or s
