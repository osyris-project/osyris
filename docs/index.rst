*********************************
Osyris - Visualization for Ramses
*********************************

Osyris is a python utility to read, manipulate and visualize simulation data
created by the astrophysical software
`Ramses <https://github.com/ramses-organisation/ramses>`_.
Osyris provides automatic handling of physical units,
loading of sub-regions inside a large simulation,
and enables the production of publication grade figures.

It was designed to be portable, lightweight and fast,
requiring minimum dependencies and resources.
It currently only works with the native ``binary`` Ramses data output format.
Osyris stores the data internally as one-dimensional arrays, and can also be used to
visualize outputs from other simulation codes, such as e.g.
`Dispatch <https://dispatch.readthedocs.io>`_.
It uses `Numpy <https://numpy.org>`_ for data
manipulation, `Pint <https://pint.readthedocs.io>`_ for physical units,
and `Matplotlib <https://matplotlib.org/stable/>`_ for visualization.

Installation
============

.. tab-set::

   .. tab-item:: pip

      .. code-block:: bash

         pip install osyris

   .. tab-item:: conda

      .. code-block:: bash

         conda install -c conda-forge osyris

   .. tab-item:: source

      .. code-block:: bash

         git clone https://github.com/osyris-project/osyris.git
         cd osyris
         python -m pip install -e .

Getting started
===============

.. grid:: 2

   .. grid-item-card:: Setup

      - :doc:`basics`
      - :doc:`configuration`

   .. grid-item-card:: Loadind data

      - :doc:`data_structures`
      - :doc:`loading_data`
      - :doc:`custom_datasets`

Plotting
========

.. grid:: 3

   .. grid-item-card:: 1D data
      :link: plotting_1d_2d.html
      :img-bottom: _images/plotting_1d_2d_5_1.png

   .. grid-item-card:: 2D data
      :link: plotting_1d_2d.html#Plotting-2D-data
      :img-bottom: _images/plotting_1d_2d_24_1.png

   .. grid-item-card:: 1D histograms
      :link: plotting_histograms.html
      :img-bottom: _images/plotting_histograms_9_1.png

.. grid:: 3

   .. grid-item-card:: 2D histograms
      :link: plotting_histograms.html#2D-histograms
      :img-bottom: _images/plotting_histograms_13_1.png

   .. grid-item-card:: Spatial maps
      :link: plotting_maps.html
      :img-bottom: _images/plotting_maps_23_1.png

   .. grid-item-card:: Thick maps
      :link: plotting_thick_maps.html
      :img-bottom: _images/plotting_thick_maps_6_1.png

.. grid:: 3

   .. grid-item-card:: Scatter plots
      :link: plotting_scatter.html
      :img-bottom: _images/plotting_scatter_11_1.png

   .. grid-item-card:: Particles
      :link: plotting_particles.html
      :img-bottom: _images/plotting_particles_5_1.png

   .. grid-item-card:: Recipes
      :link: recipes.html
      :img-bottom: _images/recipes_3_1.png

.. toctree::
    :hidden:
    :maxdepth: 2

    installation
    basics
    configuration
    data_structures
    loading_data
    custom_datasets
    plotting_1d_2d
    plotting_histograms
    plotting_maps
    plotting_thick_maps
    plotting_scatter
    plotting_particles
    recipes
    api
