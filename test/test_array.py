import osyris
import numpy as np


def test_1d_array():
    a = np.arange(100.)
    array = osyris.Array(values=a, unit="m")
    assert array.unit == osyris.units('m')
    assert len(array) == len(a)
    assert array.shape == a.shape
    assert np.allclose(array.values, a)
