

<table><tr>
  <td>
    <h1>Osyris</h1>

A python visualization utility for RAMSES astrophysical simulations data.
Osyris aims to remain portable, lightweight and fast,
to allow users to quickly explore and understand their simulation data,
as well as produce publication grade figures.
  </td>
  <td><img src="https://github.com/nvaytet/osyris/blob/main/docs/images/logo_osyris.png" width="200" /></td>
  </tr>
</table>

### Documentation ###

The documentation for `osyris` is hosted on [Readthedocs](https://osyris.readthedocs.io/en/latest/index.html).

### Installation ###

```sh
pip install osyris
```

### A short example ###

You can download the sample data [here](http://project.esss.dk/owncloud/index.php/s/biNBruU0wDOybsb/download).

Plot a 2D histogram of the cell magnetic field versus the gas density.

```python
import osyris
data = osyris.load(71, scale="au", path="data")
osyris.plot.histogram(data["density"], data["B_field"],
                      norm="log", loglog=True)
```

Create a 2D gas density slice 100 au wide through the plane normal to ``z``,
with velocity vectors overlayed as arrows, once agains using ``layers``:

```python
osyris.plane({"data": data["density"], "norm": "log"}, # layer 1
             {"data": data["velocity"], "mode": "vec"}, # layer 2
             dx=50,
             origin=center,
             direction="z")
```

### Have a problem or need a new feature? ###

Submit an issue on [Github](https://github.com/nvaytet/osyris/issues).

### Contributors ###

* Neil Vaytet (StarPlan/NBI)
* Tommaso Grassi (StarPlan/NBI)
* Matthias Gonzalez (CEA Saclay)
* Troels Haugbolle (StarPlan/NBI)
* Lucas Beeres

### Logo credit ###

[Icon vector created by frimufilms - www.freepik.com](https://www.freepik.com/free-photos-vectors/icon)
