import pathlib

import pytest

def test_import():
    """Ensure the package imports correctly"""

    import molusc

def test_path():
    """ Make sure pytest is fine with the whole package path bit
    """

    import molusc
    repo_path = pathlib.Path(molusc.__file__).resolve().parent.parent
