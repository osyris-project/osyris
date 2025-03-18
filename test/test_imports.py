# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2024 Osyris contributors (https://github.com/osyris-project/osyris)


def test_star_import_all():
    """
    Tests if the package can be imported with a star import.
    If failing, check the __init__.py files.
    """
    exec("from osyris import *", {})


def test_star_import_submodules():
    """
    Tests if the package submodules can be imported with a star import.
    If failing, check the __init__.py files.
    """
    exec("from osyris.config import *", {})
    exec("from osyris.core import *", {})
    exec("from osyris.io import *", {})
    exec("from osyris.plot import *", {})
    exec("from osyris.spatial import *", {})
    exec("from osyris.units import *", {})
