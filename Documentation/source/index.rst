MOlecular Dynamics Over Ip
=================================

What is this project for?
-------------------------

MODOI is an implementation of the Birkhoff curve shortening procedure using a technique devised by Schwetlick-Zimmer to compute local geodesics. It is currently implemented in Python using multiprocessing to distribute the computing.


How do I get set up?
----------------------------

* Install Anaconda (http://continuum.io/downloads)
* Install ASE using **pip install python-ase**
* Install pyCharm (https://www.jetbrains.com/pycharm/)
* Open pyCharm and select **Check out from Version Control**, then select **Git**
* Use Git Repository URL: https://<enter-bitbucker-username>@bitbucket.org/molecularsimulationbath/modoi.git
* Enter your BitBucket password and a copy of the project will be made locally
* The project should always be left in a state where a short simulation can be run. At the time of writing this is an all atom EMT simulation of Butane between the trans and gauche conformations. On downloaded test the software by running the **Local_Simulation.py** script.

Who do I talk to?
-----------------------

* For questions about the project please contact Daniel Sutton (sutton.c.daniel@gmail.com)

* For scientific queries please contact Johannes Zimmer (zimmer@maths.bath.ac.uk)

Components
-----------------

This software package consists of four modules:

.. toctree::
   :maxdepth: 4

   Algorithm
   Tutorial
   Implementations
   Components


Index 
==================

* :ref:`genindex`

