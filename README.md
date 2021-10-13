<table><tr><td>
<h1>Osyris</h1>
A python visualization utility for RAMSES astrophysical simulations data.
Osyris aims to remain portable, lightweight and fast,
to allow users to quickly explore and understand their simulation data,
as well as produce publication grade figures.
</td><td>
<img src="https://github.com/nvaytet/osyris/blob/main/docs/images/logo_osyris.png" width="200" />
</td></tr></table>

### Documentation ###

The documentation for `osyris` is hosted on
[Readthedocs](https://osyris.readthedocs.io).

### Installation ###

```sh
pip install osyris
```

### A short example ###

You can download the sample data
[here](https://github.com/nvaytet/osyrisdata/archive/refs/heads/main.zip).

Plot a 2D histogram of the cell magnetic field versus the gas density.

```python
import osyris
data = osyris.Dataset(8, scale="au", path="data").load()
osyris.histogram2d(data["hydro"]["density"], data["hydro"]["B_field"],
                   norm="log", loglog=True)
```
![Screenshot at 2021-10-09 21-19-12](https://user-images.githubusercontent.com/39047984/136671452-c8862fa3-fea8-421d-b3d2-630e0165fc41.png)

Create a 2D gas density slice 2000 au wide through the plane normal to ``z``,
with velocity vectors overlayed as arrows, once agains using ``layers``:

```python
ind = np.argmax(data["hydro"]["density"])
center = data["amr"]["xyz"][ind.values]
osyris.plane({"data": data["hydro"]["density"], "norm": "log"}, # layer 1
             {"data": data["hydro"]["velocity"], "mode": "vec"}, # layer 2
             dx=2000 * osyris.units("au"),
             origin=center,
             direction="z")
```
![Screenshot at 2021-10-09 21-20-01](https://user-images.githubusercontent.com/39047984/136671480-d9f16d1c-d9c9-4201-8d53-5891fae690bc.png)

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
