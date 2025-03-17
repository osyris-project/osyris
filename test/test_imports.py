# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2024 Osyris contributors (https://github.com/osyris-project/osyris)


from osyris import *  # noqa: F401, F403
from osyris.config import *  # noqa: F401, F403
from osyris.core import *  # noqa: F401, F403
from osyris.io import *  # noqa: F401, F403
from osyris.plot import *  # noqa: F401, F403
from osyris.spatial import *  # noqa: F401, F403
from osyris.units import *  # noqa: F401, F403


def test_imports_star():
    """
    Tests that the package and its submodules can be imported with a star import.
    If failing, check the __init__.py files.
    """
    pass
