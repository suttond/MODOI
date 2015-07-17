import math

import numpy as np
from scipy.linalg import orth

def mass_norm(vector, mass_matrix):
    """ This function computes the inner product of vector with mass_matrix times vector.

    Args:
      vector (numpy.array): The curve for which the length is to be computed.
      mass_matrix (numpy.array): A diagonal NumPy array containing the masses of the molecular system as computed in the SimulationClient object.

    Returns:
      float: The value of <vector, mass_matrix * vector>

    """
    return math.sqrt(np.inner(vector, mass_matrix.dot(vector)))


def orthonormal_tangent_basis(tangent, dimension):
    """ This function computes the inner product of vector with mass_matrix times vector.

    Args:
      tangent (numpy.array): The tangent direction along the local geodesic.
      dimension (int): The dimension of the problem. Computed from the atomistic simulation environment.

    Returns:
      numpy.array: An orthonomal matrix with the first column parallel to tangent.

    """

    # Set the first column of our output matrix as tangent
    mx = tangent

    # Find the first non-zero entry of the tangent vector (exists as start and endpoints are different)
    j = np.nonzero(mx)[0][0]

    # For the remaining dim - 1 columns choose unit basis vectors of the form (0,...,0,1,0,...,0) with the nonzero entry
    # not in position j.
    for i in xrange(1, dimension):
        if j != i:
            e = np.zeros(dimension)
            e[i] = 1
            mx = np.vstack((mx, e))

    mx = mx.transpose()

    # With the resulting matrix, perform the Gram-Schmidt orthonormalisation procedure on the transpose of the matrix
    # and return it.
    m, n = np.shape(mx)
    Q = np.zeros([m, n])
    R = np.zeros([n, n])
    v = np.zeros(m)

    for j in range(n):

        v[:] = mx[:,j]
        for i in range(j):
            r = np.dot(Q[:,i], mx[:,j]); R[i,j] = r
            v[:] = v[:] - r*Q[:,i]
        r = np.linalg.norm(v); R[j,j]= r
        Q[:,j] = v[:]/r

    return Q


def shifts_to_curve(start_point, end_point, shift_points, number_of_inner_points, basis_rotation_matrix,
                    tangent_direction, codimension):

    """ This function produces a curve in N-dimensional space when it is initially described as a graph over the first co-ordinate direction.

    Args:
      start_point (numpy.array): The first end point of the curve.
      end_point (numpy.array): The last end point of the curve.
      shift_points (numpy.array): Given a point x_N along the curve, the value shift_points[i] is x_N - <x_N, tangent_direction>tangent_direction. That is the orthogonal component of x_N when projected against the tangent_direction.
      number_of_inner_points (int): The number of nodes along the curve, less the end points.
      basis_rotation_matrix (numpy.array): The matrix computed as a result of the orthogonal_tangent_basis function.
      tangent_direction (numpy.array): The tangent direction as computed by the SimulationClient.
      codimension (int): The dimension of the problem minus 1. Computed from the atomistic simulation environment.

    Returns:
      list: A list of NumPy arrays describing the position of the curve.

    """

    # Compute tangent direction of line joining start and end points
    tangent = np.subtract(end_point, start_point)/(number_of_inner_points+1)

    # Initialise list to store points
    points = []

    # Generate points that are uniformly distributed along the initial line
    for i in xrange(number_of_inner_points+2):
        points.append(np.add(start_point, float(i)*tangent))

    # Shift the points as encoded in x
    for i in xrange(0, len(shift_points)/codimension):

        # Embed vector i into co_dimension + 1 dimensional space
        unrotated_shift = np.hstack((np.zeros(1), shift_points[i*codimension:(i+1)*codimension]))

        # Convert vector, by rotation, from shift from e_1 basis direction to shift from tangent direction
        shift = basis_rotation_matrix.dot(unrotated_shift)

        # Append point to list
        points[i+1] = np.add(points[i+1], shift)

    return points
