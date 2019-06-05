Faqs
=======


Installation
--------------

On **ubuntu**, you need libz-dev and python3-dev libraries which are not necessarily present by default::

  sudo apt-get install libz-dev python3-dev

Plink
~~~~~

If you have installed plink1.9 but **bioconvert** still can not use plink. It is maybe because
**bioconvert** try to call the programme by the name "plink" so you have to make a symbolic link.
First, you have to go in the repository where is plink, then use the command which: ::

    which plink1.9

go into the repository then: ::

    ln -s plink1.9 plink

after this bioconvert will be able to call plink

Libraries
-----------

Graphviz related
~~~~~~~~~~~~~~~~

When you install the requirement for dev. It is possible to have problems with pygraphviz,
try to install these libraries: ::

    sudo apt-get install libcgraph6 libgraphviz-dev

and then you can try again to install pygraphviz by two methods: ::

    pip install pygraphviz

or ::

    pip install -r requirements_dev.txt

