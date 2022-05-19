# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
"""
Define default values so that you don't have to specify them every time.
"""

parameters = {
    'path': None,
    'select': None,
    'cmap': 'viridis',
    'render_mode': 'pcolormesh',
    'sortby': {
        'part': 'identity'
    },
    'units': {},
    'constants': {}
}


def additional_variables(data):
    """
    Here are some additional variables that are to be computed every time data
    is loaded.

    It is recommended to place your variables in a `try/except` block, which
    will prevent errors if the variables are not found, for instance when
    loading data from a different simulation.
    """

    # Magnetic field
    try:
        data['hydro']['B_field'] = 0.5 * (data['hydro']['B_left'] +
                                          data['hydro']['B_right'])
    except KeyError:
        pass

    # Mass
    try:
        data['hydro']['mass'] = (data['hydro']['density'] *
                                 data['amr']['dx']**3).to('M_sun')
    except KeyError:
        pass
