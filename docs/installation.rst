***********
Instalation
***********

.. attention::

   If you are upgrading from an old version to a version >= 2.6.0, and you are getting some strange errors,
   you may need to update your ``config_osyris.py`` configuration file in ``/home/user/.osyris``.

   - If you had never touched the configuration file, it is safe to simply delete it (a new one will be created when importing ``osyris``).
   - If you had made changes to the configuration, the easiest is probably to move it to a new location/filename. Then import ``osyris`` and update the newly created file to incorporate the changes you had made previously.

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

Osyris is also available on ``conda-forge``. To install using ``conda``, do

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
