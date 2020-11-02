Optimization model of a Compressed Air Energy Storage (CAES)
===============================================

1. Installation
-----------------------------

- Get a Python distribution. See `here <http://www.anaconda.org>`_ for the latest Anaconda version and chose
  one for Python 3.x (64 or 32 bit). Run the installer afterwards.
- Download and install a `git client <https://git-scm.com/>`_ and clone this repository into a folder on your harddisk using the following command:

.. code:: bash

    git clone https://github.com/carpescientiam/CAES-optimizsation-model.git
   
A status message should confirm that all files have been downloaded. Alternatively, the files can be downloaded and extracted to a respective folder.

- Open ``Anaconda Prompt`` and type the following:

.. code:: bash

    pip install -e "C:\path\to\cloned\repository"

A message should confirm the successful installation. Please install all required packages if needed.


2. Application use
-------------------

The command lines as described below will provide some guidance how to use the model. It will be shown how you create a plant and which methodes the
created plant provides.

.. code:: python

    from caes.models import*


    # create a diabatic operating plant "caes" with the desired configuration and cavern capaicity. 
    caes = Diabatic(V_cas=310000, P_cmp=60, P_exp=321, recuperation=False)
  
  
    # obtain coefficents for linear model
    coefficients_lp1 = CAES.coefficents_linear_model()

    # change setup and re-run calculation
    CAES.P_cmp = 100
    coefficients_lp1 = CAES.coefficents_linear_model()

    # fit model for linear representation
    linear_fit = CAES.fit_cmp()
