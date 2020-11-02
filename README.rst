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

The command lines as described below will provide some guidance how to use the model. It shows how you create a plant and which methodes the
created plant provides.

To creat a diabatic operating plant "caes" without a recuperator and a compression power of 60 MW, an expansion power of 321 MW and a cavern capacity of 310000 m^3 write: 

.. code:: python

    from caes.model import*
 
    caes = Diabatic(V_cas=310000, P_cmp=60, P_exp=321, recuperation=False)
 
Methodes of the compression train:
   
.. code:: python
   
     caes.massflow_cmp(power, pressure) # returns the air massflow rate for given power and cavern presssure
     caes.temperature_cmp(power, pressure) # returns Tc1 and Tc2 for given power and cavern presssure
     caes.q_cmp(power, pressure) # returns the compression heat and the coefficient e for given power and cavern presssure
     caes.eta_exergy_cmp(power, pressure) #returns the exergetical compression efficiency

Methodes of the compression train:

.. code:: python

     caes.massflow_exp(power, pressure) 
     caes.temperature_exp(power, pressure) 
     caes.q_comb(power, pressure) # returns the required combustion heat 
     caes.eta_exergy_exp(power, pressure) 
     caes.power_heat_ratio(power, pressure)# returns the ratio: recuperation heat / expansion power
  
3. Optimization
-------------------

In order to run a operation optimization you need additionally two files as provided as "cost_example" and "prices_example". The cost_example file should carry the operation costs. The optimizsation performs based on those provided two files. To run the optimization you need to run the following command:

.. code:: python

     caes.optimize('prices_example', 'cost_example', grid_num)

The optimization method firstly executes a polynomic linear fit of the plan characteristics before the optimization itself. The name of your cost and preice files are given as a string and grid_num is a integer defining the nodes of your fit. To get the polynomic coefficients of your fit run the command:

.. code:: python

     caes.coefficents_linear_model(grid_num)
