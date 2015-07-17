Local_Simulation Script
=======================

This is a script that is included in the distribution of MODOI to run a local simulation, that is on a single machine.
It runs a single instance of SimulationServer and two instances of SimulationClient. Each instance of SimulationClient
has two instances of SimulationPotential. This particular configuration runs comfortably on a machine with 4 cores. The
script automatically handles the calculation of port numbers in this case.

.. literalinclude:: ../../Local_Simulation.py