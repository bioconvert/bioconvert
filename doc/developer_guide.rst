.. _developer_guide:

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
        _default_method = "v1"

        def __init__(self, infile, outfile):
            """
            :param str infile: information
            :param str outfile: information
            """
            super().__init__(infile, outfile)

        def _method_v1(self, *args, **kwargs):
            Conversion is made here.
            You can use self.infile and  self.outfile
            If you use an external command, you can use:
            self.execute(cmd)

        def _method_v2(self, *args, **kwargs):
            another method

You may also use this standalone to create the bioconvert_init standalone. For
instance to create the *sam* to *bam* conversion, redirect the output of the following command in
the correct file::

    $ bioconvert_init -i sam -o bam > sam2bam.py

Of course, you will need to edit the file to add the conversion itself in the
appropriate method (e.g. _method_samtools).


How to add a new method
~~~~~~~~~~~~~~~~~~~~~~~~~~

As shown above, use this coding::

    def _method_YOUuniqueMETHODname(self, *args, **kwargs):
        # from kwargs, you can use any kind of arguments.
        # threads is an example, reference, another example.
        Your code here

Then, it will be available in the class and bioconvert standalone !

The code that you will add may be of different kind:

- pure Python: just write it.
- Python code but relying on third-party library, two options:

  - if the Python library is on pypi and is simple, add it to requirements.txt
  - if the Python library requires lots of compilation, add it to requirements_tools.txt (assuming it is on bioconda).
- if the code is not on pypi or bioconda (e.g., GO code), use the self.install_tool(NAME) and add a script in ./misc/install_NAME.sh




Method decorators
~~~~~~~~~~~~~~~~~

`Decorators
<https://en.wikipedia.org/wiki/Python_syntax_and_semantics#Decorators>`_ have
been defined in ``bioconvert/core/decorators.py`` that can be used to "flag" or
"modify" conversion methods (actually, a new method is usually returned):

- ``@in_gz`` can be used to indicate that the method is able to transparenly
  handle input files that are compressed in ``.gz`` format. This is done by
  adding an ``in_gz`` attribute (set to ``True``) to the method.

- ``@compressor`` will wrap the method in code that handles input decompression
  from ``.gz`` format and output compression to ``.gz``, ``.bz2`` or ``.dsrc``.
  This automatically applies ``@in_gz``. Example:

::

    @compressor
    def _method_noncompressor(self, *args, **kwargs):
        """This method does not handle compressed input or output."""
        pass
    # This results in a method that handles compressed input and output
    # The method has an in_gz attribute (which is set to True)


- ``@out_compressor`` will wrap the method in code that handles output
  compression to ``.gz``, ``.bz2`` or ``.dsrc``. It is intended to be used on
  methods that already handle compressed input transparently, and therefore do
  not need the input decompression provided by ``@compressor``. Typically, one
  would also apply ``@in_gz`` to such methods. In that case, ``@in_gz`` should
  be applied "on top" of ``@out_compressor``. The reason is that decorators
  closest to the function are applied first, and applying another decorator on
  top of ``@in_gz`` would typically not preserve the ``in_gz`` attribute.
  Example:

::

    @in_gz
    @out_compressor
    def _method_incompressor(self, *args, **kwargs):
        """This method already handles compressed .gz input."""
        pass
    # This results in a method that handles compressed input and output
    # This method is further modified to have an in_gz attribute
    # (which is set to True)


(For more general explanations about decorators, see
https://stackoverflow.com/a/1594484/1878788.)

How to add a test and test file
-----------------------------------

Go to  ./test and add a file named ``test_fastq2fasta.py``


::

    import pytest

    from bioconvert.fastq2fasta import Fastq2Fasta
    from bioconvert import bioconvert_data
    from easydev import TempFile, md5

    @pytest.mark.parametrize("method", Fastq2Fasta.available_methods)
    def test_fastq2fasta(method):
        #your code here
        # you will need data for instance "mydata.fastq and mydata.fasta".
        # Put it in bioconvert/bioconvert/data
        # you can then use ::
        infile = bioconvert_data("mydata.fastq")
        expected_outfile = bioconvert_data("mydata.fasta")
        with TempFile(suffix=".fasta") as tempfile:
            converter = Fastq2Fasta(infile, tempfile.name)
            converter(method=method)

            # Check that the output is correct with a checksum
            assert md5(tempfile.name) == md5(expected_outfile)


Files used for testing should be added in
./bioconvert/data/testing/converter_name.
For instance test files for the
sam2paf converter should be added in
bioconvert/data/testing/sam2paf directory where you should have the test files,
a __init__.py file, a README.rst file. The latter should contain the name of the
test files and a short description.


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

you can also use the **bioconvert** standalone with -b option.


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

