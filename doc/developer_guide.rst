

.. _developer_guide:

Developer guide
===============

.. contents::

Quick start
-----------

As a developer, assuming you have a valid environment and installed Bioconvert (:ref:`install_dev`), 
go to bioconvert directory and type the bioconvert init command for the input and output formats you wish to add (here we want to convert format A to B). You may also just copy an existing file::

    cd bioconvert
    bioconvert_init -i A -o B > A2B.py

see :ref:`add_converter` section for details. Edit the file, update the method that performs the conversion by adding the relevant code (either python or external tools). Once done, please

1. add an input test file in the ./test/data directory (see :ref:`add_test`) 
2. add the relevant data test files in the ``./bioconvert/test/data/`` directory (see :ref:`add_test_file`)
3. Update the documentation as explained in :ref:`update_doc` section:

   1. add the module in doc/ref_converters.rst in the autosummary section
   2. add the A2B in the README.rst

4. add a CI action in .github/workflows named after the conversion (A2B.yml)

Note also that a converter (a Python module, e.g., fastq2fasta) may have several methods included and it is quite straightforward to add a new method (:ref:`add_method`). They can later be compared thanks to our benchmarking framework.


If this is a new formats, you may also update the glossary.rst file in the documentation.


.. _install_dev:

Installation for developers
---------------------------

To develop on `bioconvert` it is highly recommended to install `bioconvert` in a virtualenv ::

    mkdir bioconvert
    cd bioconvert
    python3.7 -m venv py37
    source py37/bin/activate

And clone the bioconvert project ::

    mkdir src
    cd src
    git clone https://github.com/bioconvert/bioconvert.git
    cd  bioconvert

We need to install some extra requirements to run the tests or build the doc so to install these requirements ::

    pip install -e . [testing]

.. warning::
    The extra requirements try to install `pygraphviz` so you need to install `graphviz` on your computer.
    If you running a distro based on debian you have to install `libcgraph6`, `libgraphviz-dev` and `graphviz` packages.

.. note::
    You may need to install extra tools to run some conversion.
    The requirements_tools.txt file list conda extra tools

.. _add_converter:

How to add a new conversion
---------------------------

Officially, **Bioconvert** supports one-to-one conversions only (from one format to
another format). See the note here below about :ref:`multi_conversion`.

Let us imagine that we want to include a new format conversion 
from :term:`FastQ` to :term:`FastA` format. 

First, you need to add a new file in the ``./bioconvert`` directory called::

    fastq2fasta.py

Please note that the name is **all in small caps** and that we concatenate the input format name, the character **2** and the output format name. Sometimes a format already includes the character 2 in its name (e.g. bz2), which may be confusing. For now, just follow the previous convention meaning duplicate the character 2 if needed (e.g., for bz2 to gz format, use bz22gz).

As for the class name, we us **all in big caps**. In the newly created file (**fastq2fasta.py**) you can (i) copy / paste the content of an existing converter (ii) use the **bioconvert_init** executable (see later), or (iii) copy / paste the following code:

.. code-block:: python
    :linenos:

    """Convert :term:`FastQ` format to :term:`FastA` formats"""
    from bioconvert import ConvBase

    __all__ = ["FASTQ2FASTA"]


    class FASTQ2FASTA(ConvBase):
        """

        """
        _default_method = "v1"

        def __init__(self, infile, outfile):
            """
            :param str infile: information
            :param str outfile: information
            """
            super().__init__(infile, outfile)

        @requires(external_library="awk")
        def _method_v1(self, *args, **kwargs):
            # Conversion is made here.
            # You can use self.infile and  self.outfile
            # If you use an external command, you can use self.execute:
            self.execute(cmd)

        @requires_nothing
        def _method_v2(self, *args, **kwargs):
            #another method
            pass


On line 1, please explain the conversion using the terms available in the :ref:`glossary`  (``./doc/glossary.rst`` file). If not available, you may edit the glossary.rst file to add a quick description of the formats.

.. warning:: If the format is not already included in **Bioconvert**, you will need to update the file core/extensions.py to add the format name and its possible extensions.

On line 2, just import the common class.

On line 7, name the class after your input and output formats; again include the
character 2 between the input and output formats. Usually, we use
big caps for the formats since most format names are acronyms. If the input or
output format exists already in **Bioconvert**, please follow the existing
conventions.

On line 13, we add the constructor.

On line 21, we add a method to perform the conversion named **_method_v1**.
Here, the prefix **_method_** is compulsary: it tells **Bioconvert** that is it a possible conversion to include in the user interface. This is also where you will add your code to perform the conversion.
The suffix name  (here **v1**) is the name of the conversion.
That way you can add as many conversion methods as you need (e.g. on line 28,
we implemented another method called **v2**).

Line 20 and line 27 show the decorator that tells **bioconvert** which external
tools are required. See :ref:`decorator` section.

Since several methods can be implemented, we need to define a default method (line 11; here **v1**).

In order to simplify the creation of new converters, you can also use the standalone **bioconvert_init**. Example::

    $ bioconvert_init -i fastq -o fasta > fastq2fasta.py

Of course, you will need to edit the file to add the conversion itself in the
appropriate method (e.g. _method_v1).



If you need to include extra arguments, such as a reference file, you may add extra argument, although this is not yet part of the official **Bioconvert** API. See for instance :class:`~bioconvert.sam2cram.SAM2CRAM` converter.



.. _multi_conversion:

One-to-many and many-to-one conversions
---------------------------------------

The one-to-many and many-to-one conversions are now implemented in
**Bioconvert**. We have only 2 instances so far namely class:`bioconvert.fastq2fasta_qual`
and  class:`bioconvert.fasta_qual2fastq`. We have no instances of many-to-many
so far. The underscore character purpose is to indicate a **and** connection. So
you need QUAL *and* FASTA to create a FASTQ file.

For developers, we ask the input or output formats to be sorted alphabetically
to ease the user experience.


.. _add_method:

How to add a new method to an existing converter
------------------------------------------------

As shown above, use this code and add it to the relevant file in ``./bioconvert``
directory::

    def _method_UniqueName(self, *args, **kwargs):
        # from kwargs, you can use any kind of arguments.
        # threads is an example, reference, another example.
        # Your code here below
        pass

Then, it will be available in the class and **bioconvert** 
automatically; the **bioconvert** executable should show the name of your new method in the help message.

In order to add your new method, you can add:

* Pure Python code
* Python code that relies on third-party library. If so, you may use:
  
    * Python libraries available on pypi. Pleaes add the library name to the
      requirements.txt
    * if the Python library requires lots of compilation and is available
      on bioconda, you may add the library name to the requirements_tools.txt
      instead.
      
* Third party tools available on **bioconda** (e.g., squizz, seqtk, etc)
  that you can add to the requirements_tools.txt
* Perl and GO code are also accepted. If so, use the self.install_tool(NAME)
  and add a script in ``./misc/install_NAME.sh``


.. _decorator:

Decorators
----------

`Decorators
<https://en.wikipedia.org/wiki/Python_syntax_and_semantics#Decorators>`_ have
been defined in ``bioconvert/core/decorators.py`` that can be used to "flag" or
"modify" conversion methods:

- ``@in_gz`` can be used to indicate that the method is able to transparently
  handle input files that are compressed in ``.gz`` format. This is done by
  adding an ``in_gz`` attribute (set to ``True``) to the method.

- ``@compressor`` will wrap the method in code that handles input decompression
  from ``.gz`` format and output compression to ``.gz``, ``.bz2`` or ``.dsrc``.
  This automatically applies ``@in_gz``.

  Example:

::

    @compressor
    def _method_noncompressor(self, *args, **kwargs):
        """This method does not handle compressed input or output by itself."""
        pass
    # The decorator transforms the method that now handles compressed 
    # input and output; the method has an in_gz attribute (which is set to True)


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


Another **bioconvert** decorator is called **requires**. 

It should be used to annotate a method with the type of tools it needs to work.

It is important to decorate all methods with the **requires** decorator so that user
interface can tell what tools are properly installed or not. You can use 4
arguments as explained in :mod:`bioconvert.core.decorators`:


.. code-block:: python
    :linenos:

    @requires_nothing
    def _method_python(self, *args, **kwargs):
        # a pure Python code does not require extra libraries
        with open(self.outfile, "w") as fasta, open(self.infile, "r") as fastq:
             for (name, seq, _) in FASTQ2FASTA.readfq(fastq):
                 fasta.write(">{}\n{}\n".format(name, seq))

     @requires(python_library="mappy")
     def _method_mappy(self, *args, **kwargs):
         with open(self.outfile, "w") as fasta:
             for (name, seq, _) in fastx_read(self.infile):
                 fasta.write(">{}\n{}\n".format(name, seq))

     @requires("awk")
     def _method_awk(self, *args, **kwargs):
         # Note1: since we use .format, we need to escape the { and } characters
         # Note2: the \n need to be escaped for Popen to work
         awkcmd = """awk '{{printf(">%s\\n",substr($0,2));}}' """
         cmd = "{} {} > {}".format(awkcmd, self.infile, self.outfile)
         self.execute(cmd)


On line 1, we decorate the method with the :func:`~bioconvert.core.decorators.requires_nothing` decorator because the method is implemented in Pure Python.

One line 8, we decorate the method with the :func:`~bioconvert.core.decorators.requires` decorator to inform **bioconvert** that the method relies on the external Python library called mappy. 


One line 14, we decorate the method with the :func:`~bioconvert.core.decorators.requires` decorator to inform **bioconvert** that the method relies on an external tool called awk. In theory, you should write::

    @requires(external_library="awk")

but ``external_library`` is the first optional argument so it can be omitted. If several libraries are required, you can use::

    @requires(external_libraries=["awk", ""])

or::

    @requires(python_libraries=["scipy", "pandas"])


.. note:: For more general explanations about decorators, see https://stackoverflow.com/a/1594484/1878788.


.. _add_test:

How to add a test
-----------------

Following the example from above (fastq2fasta), we need to add a test file. To
do so, go to the ``./test`` directory and add a file named ``test_fastq2fasta.py``.

.. code-block:: python
    :linenos:


    import pytest

    from bioconvert.fastq2fasta import FASTQ2FASTA
    from bioconvert import bioconvert_data
    from easydev import TempFile, md5

    from . import test_dir

    @pytest.mark.parametrize("method", FASTQ2FASTA.available_methods)
    def test_fastq2fasta(method):
        # your code here
        # you will need data for instance "mydata.fastq and mydata.fasta".
        # Put it in bioconvert/bioconvert/data
        # you can then use ::
        infile = f"{test_dir}/data/fastq/test_mydata.fastq"
        expected_outfile = f"{test_dir}/data/fasta/test_mydata.fasta"
        with TempFile(suffix=".fasta") as tempfile:
            converter = FASTQ2FASTA(infile, tempfile.name)
            converter(method=method)

            # Check that the output is correct with a checksum
            assert md5(tempfile.name) == md5(expected_outfile)


In **Bioconvert**, we use **pytest** as our test framework. In principle, we 
need one test function per method found in the converter. Here
on line 7 we serialize the tests by looping through the methods available in the
converter using the pytest.mark.parametrize function. That way, the test 
file remains short and do not need to be duplicated.



.. _add_test_file:

How to add a test file
----------------------

Files used for testing should be added in ``./bioconvert/test/data/ext/converter_name.ext``.

How to locally run the tests
----------------------------

Go to the source directory of **Bioconvert**. 

If not already done, install all packages required for testing::

    cd bioconvert
    pip3 install .[testing]

Then, run the tests using::

    pytest test/ -v

Or, to run a specific test file, for example for your new converter fastq2fasta::

    pytest test/test_fastq2fasta.py -v

or ::

    pytest -v -k test_fastq2fasta


How to benchmark your new method vs others
------------------------------------------

::

    from bioconvert import Benchmark
    from bioconvert.fastq2fasta import FASTQ2FASTA
    converter = FASTQ2FASTA(infile, outfile)
    b = Benchmark(converter)
    b.plot()

you can also use the **bioconvert** standalone with -b option.


.. _update_doc:

How to add you new converter to the main documentation ?
--------------------------------------------------------

Edit the doc/ref_converters.rst and add this code (replacing A2B by your conversion)::


    .. automodule:: bioconvert.A2B
        :members:
        :synopsis:
        :private-members:

and update the autosummary section::

    .. autosummary::

        bioconvert.A2B


pep8 and conventions
--------------------

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


Since v0.5.2, we apply black on the different Python modules.

Requirements files
------------------

- requirements.txt : should contain the packages to be retrieved from Pypi only.
  Those are downloaded and installed (if missing) when using
  **python setup.py install**
- environment_rtd.yml : do not touch. Simple file for readthedocs
- readthedocs.yml : all conda and pip dependencies to run the example and build
  the doc
- environment.yml is a conda list of all dependencies

How to update bioconvert on bioconda
------------------------------------

Fork bioconda-recipes github repository and clone locally. Follow instructions on
https://bioconda.github.io/contributing.html

In a nutshell, install bioconda-utils::

    git clone YOURFORKED_REPOSITORY
    cd bioconda-recipes

edit bioconvert recipes and update its contents. If a new version pypi exists, you need to change the md5sum in ``recipes/bioconvert/meta.yaml``.


check the recipes::

    bioconda-utils build  recipes/ config.yml --packages bioconvert

Finally, commit and created a PR::

    #git push -u origin my-recipe
    git commit .
    git push


Sphinx Documentation
--------------------

In order to update the documentation, go the *./doc* directory and update any of
the .rst file. Then, for Linux users, just type::

    make html

Regarding the :ref:`formats` page, we provide simple ontology with 3 entries:
Type, Format and Status. Please choose one of the following values:

- Type: sequence, assembly, alignement, other, index, variant, database,  compression
- Format: binary, human-readable
- Status: deprecated, included, not included


Docker
------
In order to create the docker file, use this command::

    docker build .

The Dockerfile found next to setup.py is self-content and has been tested for v0.5.2 ; it uses the spec-file.txt that was generated in a conda environment using::

    conda list --explicit

