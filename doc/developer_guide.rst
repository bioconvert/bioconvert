Developer guide
=================


How to add a new converter ?
-----------------------------------

For now, converters are simple conversion from one format to another one.
There is no third-party file. For instance, if you need a reference file, this
is not part of the API for the moment.

Now, let us take a simple example such as a fastq to fasta conversion.

First, you need to create a module (Python file). We use the convention::

    input2output.py

all in small caps ! In this file, copy and paste this example::


    """Convert :term:`FastQ` format to :term:`Fasta` formats"""
    from bioconvert import ConvBase

    __all__ = ["Fastq2Fasta"]


    class Fastq2Fasta(ConvBase):
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
            self._default_method = "v1"

        def _method_v1(self, *args, **kwargs):
            Conversion is made here.
            You can use self.infile and  self.outfile
            If you use an external command, you can use:
            self.execute(cmd)

        def _method_v2(self, *args, **kwargs):
            another method

You may also use this standalone to create the bioconvert_init standalone. For
instance to create the bz2 to gz convertm redirect the output of this command in
the correct file::

    $ bioconvert_init -i bz2 -o gz > bz22gz.py

Of course, you will need to edit the file to add the conversion itself in the
appropriate method (e.g. _method_gz).


How to add a test
-----------------------

Go to  ./test and add a file named **test_fastq2fasta.py**


::

    import pytest
    from bioconvert.fastq2fasta import Fastq2Fasta

    def test_fastq2fasta():
        #your code here
        # you will need data for instance "mydata".
        # Put it in bioconvert/bioconvert/data
        # you can then use ::
        from bioconvert import bioconvert_data
        bioconvert_data("mydata")


How to locally run the tests
----------------------------

Go to root directory. If not already done, install all packages listed in ``requirements_dev.txt``.
You can do so by running::

    pip3 install -r requirements_dev.txt

Then, run the tests using::

    pytest test/ -v

Or, to run a specific test file, for example for your new convertor fastq2fasta::

    pytest test/test_fastq2fasta.py -v



How to benchmark your new method vs others
--------------------------------------------------

::

    from bioconvert import Benchmark
    from bioconvert.fastq2fasta import Fastq2Fasta
    converter = Fastq2Fasta(infile, outfile)
    b = Benchmark(converter)
    b.plot()



How to add you new converter to the main documentation ?
-----------------------------------------------------------

Edit the doc/references.rst and add those lines ::

    .. automodule:: bioconverter.fastq2fasta
        :members:
        :synopsis:


pep8 and conventions
-------------------------

In order to write your Python code, use PEP8 convention as much as possible.
Follow the conventions used in the code. For instance,

::

    class A():
        """Some documentation"""

        def __init__(self):
            """some doc"""
            pass

        def another_method(self):
            """some doc"""
            c = 1 + 2


    class B():
        """Another class"""

        def __init__(self, *args, **kwargs):
            """some doc"""
            pass


     def AFunction(x):
        """some doc"""
        return x


- 2 blank lines between  classes and functions
- 1 blank lines between methods
- spaces around operators (e.g. =, +)
- Try to have 80 characters max on one line
- Add documentation in triple quotes


To check PEP8 compliance of a python source code file, you can run ``flake8`` on it.
For instance::

    $ flake8 bioconvert/fastq2fasta.py

Requirements files
------------------------

- requirements.txt : should contain the packages to be retrieved from Pypi only.
  Those are downloaded and installed (if missing) when using
  **python setup.py install**
- environment_rtd.yml : do not touch. Simple file for readthedocs
- readthedocs.yml : all conda and pip dependencies to run the example and build
  the doc
- requirements_dev.txt : packages required for testing or building the doc (not
  required to run the bioconvert package
- requirements_tools.txt : all conda dependencies

