.. osyris documentation master file, created by
   sphinx-quickstart on Mon Jul 15 16:29:19 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

********************
Osyris documentation
********************

Osyris is a python visualization utility for RAMSES data. Its purpose is to
plot quick diagnostics while a simulation is running, and also produce
publication grade figures.

Osyris was developed to provide a light-weight method to read
`RAMSES <https://bitbucket.org/rteyssie/ramses>`_ simulation outputs using
``python``.
It currently only works with the native ``binary`` data output format.
It uses ``numpy`` for data manipulation and ``matplotlib`` for visualization.

Contents
========

.. toctree::
   :maxdepth: 2
   :caption: Getting started:

   getting-started/installation
   getting-started/load-ramses-data
   getting-started/2d-histogram
   getting-started/2d-slice

.. toctree::
   :maxdepth: 2
   :caption: Demos:

   demos

.. toctree::
   :maxdepth: 2
   :caption: API:

.. toctree::
   :maxdepth: 2
   :caption: Support:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
