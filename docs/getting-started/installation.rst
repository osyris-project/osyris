.. _installation:

Installation
============

Using ``pip``
-------------

.. code-block:: sh

   $ pip install osyris


Clone the source
----------------

**Requirements**: you will need `matplotlib` installed on your system.
Clone the ``osyris`` repository from [Github](https://github.com/nvaytet/osyris)
into your chosen directory. For this tutorial, the directory will be located
at ``/home/user/software``:

.. code-block:: sh

   $ cd /home/user/software
   $ git clone https://github.com/nvaytet/osyris.git


Then just add the path to the ``src/osyris`` directory to you ``PYTHONPATH``:

.. code-block:: sh

   $ export PYTHONPATH=$PYTHONPATH:/home/user/software/osyris/src/osyris

To avoid having to do this every time you open a new terminal, you can of course
add this line to your ``.bashrc`` file.
