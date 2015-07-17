The Algorithm
=============

The purpose of this code is to provide a consistent numerical method to find the minimum of the functional

.. math::

   L[u] = \int_0^1\sqrt{2(E-V(u(\tau)))}\|u'(\tau)\| d\tau

over the class of Lipschitz curves joining :math:`q_a` to :math:`q_b`  in :math:`\mathbb R^d`. Here we will assume that
:math:`V \in C^2(\mathbb R^d)` and :math:`E > \| V \|_{\infty}`. The Maupertuis principle states that by computing this
minimum, we consequently compute solutions to

.. math::

   q''(t) = -\nabla V(q(t))

subject to :math:`q(0) = \xi_1` and :math:`q(T) = \xi_2` for some :math:`T` to be determined as a function of :math:`E`.
The MODOI software package uses effective medium theory to compute the value of V based on the atomistic description of
the system.

For full details consult Chapter 4 of Microscopic Hamiltonian Systems and their Effective Description :download:`pdf <DC_Sutton_Thesis.pdf>`.