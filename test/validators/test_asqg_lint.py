import pytest
from bioconvert.validator.asqg_lint import ASQGLint
from bioconvert import bioconvert_data



def test_asqg_lint():


    for filename in ['test_wrong1.asqg']:
        data = bioconvert_data(filename, "validator")
        try:
            ASQGLint(data).validate()
            assert False
        except:
            assert True

    for filename in ['test_example1.asqg', "test_empty_lines.asqg"]:
        data = bioconvert_data(filename, "validator")
        ASQGLint(data).validate()



