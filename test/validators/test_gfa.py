import pytest
from bioconvert.validator.gfa_lint import GFALint
from bioconvert import bioconvert_data



def test_gfa_lint():


    for filename in ['test_gfa2fasta.gfa']:
        data = bioconvert_data(filename, "validator")
        GFALint(data).validate()



