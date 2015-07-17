import numpy as np
import math

def norm(x, matrix):
    """ Computes the value of sqrt(<x, matrix*x>).

    Args:
      x (numpy.array) :
          A vector, stored as a NumPy array, to compute the norm for.
      matrix (numpy.array) :
          A matrix, stored as a NumPy array, used in the computation of <x, matrix*x>.


    Returns:
      float :
          The value of sqrt(<x, matrix*x>).

    """

    return math.sqrt(np.inner(x, matrix.dot(x)))


def norm_gradient(x, matrix):
    """ Computes the gradient of sqrt(<x, matrix*x>).

    Args:
      x (numpy.array) :
          A vector, stored as a NumPy array, to compute the norm for.
      matrix (numpy.array) :
          A matrix, stored as a NumPy array, used in the computation of <x, matrix*x>.


    Returns:
      numpy.array :
          The gradient of sqrt(<x, matrix*x>).

    """

    a = (matrix + matrix.transpose())

    return a.dot(x) / (2 * norm(x, matrix))


def Length(curve, metric, number_of_inner_nodes, mass_matrix):
    """ This function computes the length of a curve object in the isotropic Riemannian length functional with metric
    coefficient metric.

    Args:
      curve (Curve): The curve for which the length is to be computed.
      metric: A list of float values representing the values of the metric along the curve. We have that metric[i] = a(curve[i]) where a is the metric coefficient.
      number_of_inner_nodes (int): The number of nodes in the curve object, less the end points.
      mass_vector (numpy.array): A NumPy array containing the masses of the molecular system as computed in the SimulationClient object.

    Returns:
      float: The length of the curve with metric values in metric.

    """

    # Convert the shifts x into points in the full dimensional space
    points = curve

    # Pre-compute the metric values to minimise repeated metric evaluations
    a = metric

    # Compute quantities used to determine the length and gradient
    n = np.subtract(points[1], points[0])
    b = norm(n, mass_matrix)
    u = (a[1][0]+a[0][0])

    # Initialise the length with the trapezoidal approximation of the first line segments length
    l = u * b

    for i in xrange(1, len(points)-1):

        # Compute the quantities needed for the next trapezoidal rule approximation.
        n = np.subtract(points[i+1], points[i])
        v = (a[i+1][0]+a[i][0])

        # Add length of line segment to total length
        l += v * norm(n, mass_matrix)

    return 0.5 * l


def GradLength(curve, metric, number_of_inner_nodes, mass_matrix, basis_rotation_matrix):
    """ This function computes the gradient of the length of a curve object in the isotropic Riemannian length
    functional with metric coefficient metric.

    Args:
      curve (Curve): The curve for which the length is to be computed.
      metric: A list of float values representing the values of the metric along the curve. We have that metric[i] = a(curve[i]) where a is the metric coefficient.
      number_of_inner_nodes (int): The number of nodes in the curve object, less the end points.
      mass_vector (numpy.array): A NumPy array containing the masses of the molecular system as computed in the SimulationClient object.

    Returns:
      numpy.array: The gradient of the length functional on the curve with metric values in metric.

    """

    # Convert the shifts x into points in the full dimensional space
    points = curve

    # Pre-compute the metric values to minimise repeated metric evaluations
    a = metric

    # Compute quantities used to determine the length and gradient
    n = np.subtract(points[1], points[0])
    b = norm(n, mass_matrix)
    c = norm_gradient(n, mass_matrix)
    u = (a[1][0]+a[0][0])

    # Initialise a list to store the gradient
    g = []

    for i in xrange(1, len(points)-1):

        # Compute the quantities needed for the next trapezoidal rule approximation.
        n = np.subtract(points[i+1], points[i])
        d = norm(n, mass_matrix)
        e = norm_gradient(n, mass_matrix)
        v = (a[i+1][0]+a[i][0])

        # Compute next gradient component and update gradient
        g.append(basis_rotation_matrix.transpose().dot(a[i][1] * (b + d) + u * c - v * e)[1:])

        # Pass back calculated values for efficiency
        b = d
        c = e
        u = v

    return 0.5 * np.asarray(g).flatten()