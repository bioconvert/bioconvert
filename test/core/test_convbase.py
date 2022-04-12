from bioconvert import ConvBase
from bioconvert.bam2cov import BAM2COV
from easydev import TempFile

from .. import test_dir


def test_convbase():

    # General tests
    infile = f"{test_dir}/data/bam/test_measles.sorted.bam"
    with TempFile(suffix=".cov") as outfile:
        c = BAM2COV(infile, outfile.name)
        c()
        try:
            c(method_name="wrong")
            assert False
        except ValueError:
            assert True

        c.shell("ls")
        c.execute("ls", ignore_errors=True, verbose=True)
        c.execute("ls", ignore_errors=True, verbose=True, shell=True)
        c._execute("ls")

    # Wrong name
    try:

        class TEST(ConvBase):
            input_ext = ".fa"
            output_ext = ".fq"

            def __call__(self):
                pass

        assert False
    except:
        assert True

    # add dot
    class in2out(ConvBase):
        input_ext = "in"
        output_ext = "out"

        def __call__(self):
            pass

    # wrong input extension (int)
    try:

        class int2out(ConvBase):
            input_ext = [1]
            output_ext = ".out"

            def __call__(self):
                pass

        assert False
    except:
        assert True

    # add dot  mix case
    class in2out(ConvBase):
        input_ext = ["in", ".in2"]
        output_ext = "out"

        def __call__(self):
            pass

    try:

        class in2out(ConvBase):
            input_ext = 1
            output_ext = 2

            def __call__(self):
                pass

        assert False
    except:
        assert True

    class in2out(ConvBase):
        input_ext = [".fa"]
        output_ext = [".fq"]

        def __call__(self):
            self.execute("ls")

    this = in2out(infile, "test.fq")
    assert this.name == "in2out"
    this()
