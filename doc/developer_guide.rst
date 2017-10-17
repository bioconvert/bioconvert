Developer guide
=================


How to add a new converter ?
-----------------------------------

What is accepted as a converter a simple conversion from one foramt to another
one without third-party file. For instance, if you need a reference file, this
is not part of the API for the moment. 

Now, let us take a simple example such as a fastq to fasta conversion. 


First, you need to create a module (Python file). We use the convention::


    input2output.py

all in small caps if possible.

In this file, copy and paste this example::


    """Convert :term:`FastQ` format to :term:`Fasta` formats"""
    from bioconvert import ConvBase

    __all__ = ["Fastq2FastQ"]


    class Fastq2FastQ(ConvBase):
        """

        """
        input_ext = ['.fastq']
        output_ext = '.fasta'

        def __init__(self, infile, outfile):
            """
            :param str infile: information
            :param str outfile: information
            """
            super().__init__(infile, outfile)

        def __call__(self, *args, **kwargs):
            Conversion is made here. 
            You can use self.infile and  self.outfile
            If you use an external command, you can use:

            self.execute(cmd)


