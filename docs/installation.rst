************
Installation
************

Using pip
=========

To install using ``pip``, simply do

.. code-block:: bash

   pip install osyris

If you have already installed ``osyris`` in the past, but would like to upgrade, use

.. code-block:: bash

   pip install osyris --upgrade

Using conda
===========

Osyris is also available on ``conda-forge``. To install using ``conda``, use

.. code-block:: bash

   conda install -c conda-forge osyris

From source
===========

To install from source, clone the ``osyris`` repository from `Github <https://github.com/osyris-project/osyris>`_
into your chosen directory.
For this tutorial, the directory will be located at ``/home/user/software``.
Then just do a local editable ``pip`` install:

.. code-block:: bash

   cd /home/user/software
   git clone https://github.com/osyris-project/osyris.git
   cd osyris
   python -m pip install -e .

Requirements
============

Osyris requires the following packages:

- `Numpy <https://numpy.org>`_
- `Matplotlib <https://matplotlib.org>`_
- `Numba <https://numba.pydata.org>`_
- `Pint <https://pint.readthedocs.io>`_
