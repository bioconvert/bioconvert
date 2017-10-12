BioKit
==========

Bioinformatics tools in Python




.. image:: https://badge.fury.io/py/biokit.svg
    :target: https://pypi.python.org/pypi/biokit

.. image:: https://secure.travis-ci.org/biokit/biokit.png
    :target: http://travis-ci.org/biokit/biokit

.. image:: https://coveralls.io/repos/cokelaer/biokit/badge.png?branch=master 
   :target: https://coveralls.io/r/cokelaer/biokit?branch=master 



:note: BioKit is tested with Travis for the following Python version: 2.7.9
       3.4.2 and 3.5.0

:contributions: Please join https://github.com/biokit/biokit and share your notebooks https://github.com/biokit/biobooks/
:issues: Please use https://github.com/biokit/biokit/issues


.. image:: http://pythonhosted.org/biokit/_images/biokit.gif
    :target: http://pythonhosted.org/biokit/_images/biokit.gif

Contents
===============

BioKit is a set of tools gathered from several other Python packages. The goal
is to gather tools that should be useful to develop computational biology
software. Biokit contains a few plotting tools (viz module), some statistical
analysis (mixture model), some tools to access to Taxon and GO identifier, some basic tools to manipulate sequences and so on. It is linked to BioServices package to provide access to biological resources. Lots of biological software are developed in R. We have also added a module to ease the installation and usage of R tools within BioKit.

Installation
==============

::

    pip install biokit


Note about testing
====================

From travis, coverage is about 50% at the moment, which is low because some tests are ignored. Tests ignored are
those that are slow or required R dependencies. To be ignored, we filled the setup.cfg with an option called **attr**. 
IF you comment that attribute in the **setup.cfg** and run ::

    python setup.py nosetests
    
You should reach a higher coverage (about 70%)    
