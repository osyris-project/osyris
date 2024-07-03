[![Documentation Status](https://readthedocs.org/projects/osyris/badge/?version=latest)](https://osyris.readthedocs.io/en/stable/?badge=latest)
[![Join the chat at https://app.gitter.im/#/room/#osyris-project_community:gitter.im](https://badges.gitter.im/Join%20Chat.svg)](https://app.gitter.im/#/room/#osyris-project_community:gitter.im)

# Osyris

A python visualization utility for RAMSES astrophysical simulations data.
Osyris aims to remain portable, lightweight and fast,
to allow users to quickly explore and understand their simulation data,
as well as produce publication grade figures.

## Documentation

The documentation for `osyris` can be found at https://osyris.readthedocs.io.

## Installation

```sh
pip install osyris
```

## A short example

You can download the sample data
[here](https://github.com/osyris-project/osyrisdata/archive/refs/heads/main.zip).

Plot a 2D histogram of the cell magnetic field versus the gas density.

```python
import numpy as np
import osyris

data = osyris.RamsesDataset(8, path="data").load()
osyris.hist2d(data["mesh"]["density"], data["mesh"]["B_field"],
              norm="log", loglog=True)
```
![hist2d](https://osyris.readthedocs.io/en/stable/_images/plotting_histograms_13_1.png)

Create a 2D gas density map 2000 au wide through the plane normal to ``z``,
with velocity vectors overlayed as arrows, once again using ``layers``:

```python
ind = np.argmax(data["mesh"]["density"])
center = data["mesh"]["position"][ind]
osyris.map(
    data["mesh"].layer("density", norm="log"),
    data["mesh"].layer("velocity", mode="vec"),
    dx=2000 * osyris.units("au"),
    origin=center,
    direction="z",
)
```
![map2d](https://osyris.readthedocs.io/en/stable/_images/plotting_maps_23_1.png)

## Have a problem or need a new feature?

- Bug reports or feature requests should be submitted by opening an [issue](https://github.com/osyris-project/osyris/issues)
- For general discussions or questions about how to do something with `osyris`, start a new [discussion](https://github.com/osyris-project/osyris/discussions)

## Logo credit

[Icon vector created by frimufilms - www.freepik.com](https://www.freepik.com/free-photos-vectors/icon)
