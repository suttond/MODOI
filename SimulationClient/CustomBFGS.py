import numpy as np
from scipy.optimize import minpack2

import LinearAlgebra as la
from Geometric import Length, GradLength
from MetricValues import get_metric


def find_geodesic_midpoint(start_point, end_point, number_of_inner_points, basis_rotation_matrix,
                           tangent_direction, codimension, metric_server_addresses, mass_matrix, authkey,
                           gtol=1e-5):
    """ This function computes the local geodesic curve joining start_point to end_point using a modified BFGS method.
    The modification arises from taking the implementation of BFGS and re-writing it to minimise the number
    of times the metric function is called.

    Args:
      start_point (numpy.array) :
          The first end point of the curve.
      end_point (numpy.array) :
          The last end point of the curve.
      number_of_inner_points (int) :
          The number of nodes along the curve, less the end points.
      basis_rotation_matrix (numpy.array) :
          The matrix computed as a result of the orthogonal_tangent_basis function.
      tangent_direction (numpy.array) :
          The tangent direction as computed by the SimulationClient.
      codimension (int) :
          The dimension of the problem minus 1. Computed from the atomistic simulation environment.
      metric_server_addresses :
          A list of tuples of the form (str, int) containing the hostnames and port numbers of the SimulationPotential
          instances.
      mass_matrix (numpy.array) :
          A diagonal NumPy array containing the masses of the molecular system as computed in the SimulationClient
          object.
      authkey (str) :
          The password used in order to communicate with the SimulationPotential instances.
      gtol (optional float) :
          The tolerance threshold for the BGFS method.

    Returns:
      numpy.array: The midpoint along the local geodesic curve.

    """

    # Determine the number of variable the BFGS method will be applied to
    number_of_variables = number_of_inner_points * codimension

    # Produce an initial guess for the minimum, in this case it will be that the straight line segment joining
    # start_point to end_point is the initial guess
    x0 = np.zeros(number_of_variables)

    # Allocate memory for the kth iterate
    xk = x0

    # Convert the description of the curve as shifts in the orthonormal hyperspace along the initial line to points in
    # the full space. See LinearAlgebra.shifts_to_curve for more details.
    curve = la.shifts_to_curve(start_point, end_point, xk, number_of_inner_points,
                            basis_rotation_matrix, tangent_direction, codimension)

    # Get the initial metric values along the starting curve.
    metric = get_metric(curve, number_of_inner_points, metric_server_addresses, authkey)

    # If the SimulationPotential couldn't be contacted then return None to close the SimulationClient
    if metric is None:
        return None

    # Obtain the initial gradient of the length functional along the curve
    gfk = GradLength(curve, metric, number_of_inner_points, mass_matrix, basis_rotation_matrix)

    # Create an identity matrix object
    I = np.eye(number_of_variables, dtype=int)

    # Initialise the memory to store the approximate Hessian matrix
    Hk = I

    # Compute the norm of the gradient in the L^{\infty} norm
    gnorm = np.amax(np.abs(gfk))

    # The main body of the BFGS calculation:
    # Repeat the method until the norm of the gradient of the length is sufficiently small.
    while gnorm > gtol:

        alpha1 = 1.0
        pk = -np.dot(Hk, gfk)

        phi0 = Length(curve, metric, number_of_inner_points, mass_matrix)
        phi1 = phi0
        derphi0 = np.dot(gfk, pk)
        derphi1 = derphi0

        isave = np.zeros((2,), np.intc)
        dsave = np.zeros((13,), float)
        task = b'START'

        # Perform the linesearch
        for i in xrange(30):
            stp, phi1, derphi1, task = minpack2.dcsrch(alpha1, phi1, derphi1, 1e-4, 0.9, 1e-14, task, 1e-8, 50,
                                                       isave, dsave)
            if task[:2] == b'FG':
                alpha1 = stp
                # Convert the description of the curve as shifts in the orthonormal hyperspace along the initial line
                # to points in the full space. See LinearAlgebra.shifts_to_curve for more details.
                curve = la.shifts_to_curve(start_point, end_point, xk + stp*pk, number_of_inner_points,
                                        basis_rotation_matrix, tangent_direction, codimension)

                # Get the initial metric values along the current trial.
                metric = get_metric(curve, number_of_inner_points, metric_server_addresses, authkey)

                # If the SimulationPotential couldn't be reached then return None to close SimulationClient
                if metric is None:
                    return None

                phi1 = Length(curve, metric, number_of_inner_points, mass_matrix)
                gfkp1 = GradLength(curve, metric, number_of_inner_points, mass_matrix, basis_rotation_matrix)
                derphi1 = np.dot(gfkp1, pk)
            else:
                break
        else:
            break

        if task[:5] == b'ERROR' or task[:4] == b'WARN':
            break

        alpha_k = stp
        xkp1 = xk + alpha_k * pk
        sk = xkp1 - xk
        xk = xkp1
        yk = gfkp1 - gfk
        gfk = gfkp1
        gnorm = np.amax(np.abs(gfk))
        if gnorm <= gtol:
            break

        rhok = 1.0 / (np.dot(yk, sk))
        if np.isinf(rhok): rhok = 1000.0  # this is patch for numpy

        Hk = np.dot(I - sk[:, np.newaxis] * yk[np.newaxis, :] *
                    rhok, np.dot(Hk, I - yk[:, np.newaxis] * sk[np.newaxis, :] * rhok)) + (rhok * sk[:, np.newaxis]
                                                                                           * sk[np.newaxis, :])

    # Return the midpoint
    return curve[(number_of_inner_points + 1) / 2]