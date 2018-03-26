import pytest
import os

skiptravis = pytest.mark.skipif(
    "TRAVIS_PYTHON_VERSION" in os.environ and
    os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"),
    reason="On travis",
)
