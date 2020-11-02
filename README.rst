Models for compressed air energy storage (CAES)
===============================================

1.Installation
-----------------------------

- Get a Python distribution. See `here <http://www.anaconda.org>`_ for the latest Anaconda version and chose
  one for Python 3.x (64 or 32 bit). Run the installer afterwards.
- Download and install a `git client <https://git-scm.com/>`_ and clone this repository into a folder on your harddisk using the following command:

.. code:: bash

    git clone https://github.com/carpescientiam/CAES-optimizsation-model.git

    A status message should confirm that all files have been downloaded.
    Alternatively, the files can be downloaded and extracted to a respective folder.

- Open ``Anaconda Prompt`` and type the following:

.. code:: bash

    pip install -e "C:\path\to\cloned\repository"

A message should confirm the successful installation.

- Install required packages
Run the application
-------------------

Run the command line application as described below.

.. code:: python

    from caes.models import*


    # create setup and calculate coefficients (default is with recuperation)
    CAES = DiabaticCAES(V_cas=310000, P_cmp=60, P_exp=321, recuperation=False)

    # obtain coefficents for linear model
    coefficients_lp1 = CAES.coefficents_linear_model()

    # change setup and re-run calculation
    CAES.P_cmp = 100
    coefficients_lp1 = CAES.coefficents_linear_model()

    # fit model for linear representation
    linear_fit = CAES.fit_cmp()
