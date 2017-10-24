from bioconvert import ConvBase
from bioconvert.bam2bed import BAM2BED
from bioconvert import bioconvert_data
from easydev import TempFile



def test_convbase():
    infile = bioconvert_data("test_measles.fa")

    with TempFile(suffix=".bed") as outfile:
        BAM2BED(infile, outfile.name)

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
    this = in2out("test.fa", "test.fq")
    assert this.name== "in2out"
    this()


