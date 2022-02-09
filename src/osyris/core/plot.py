# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/nvaytet/osyris)


class Plot:
    def __init__(self, x=None, y=None, layers=None, fig=None, ax=None, filename=None):
        self.x = x
        self.y = y
        self.layers = layers
        self.fig = fig
        self.ax = ax
        if filename is not None:
            self.fig.savefig(filename, bbox_inches='tight')
