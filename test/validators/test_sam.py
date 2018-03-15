import pytest
from bioconvert.validator.sam_lint import SAMLint
from bioconvert import bioconvert_data



def test_gfa_lint():

    filename = bioconvert_data("test_measles.sam")
    SAMLint(filename).validate()



