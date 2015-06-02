import numpy as np

import LinearAlgebra as la


def Length(curve, metric, number_of_inner_nodes, mass_vector):
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

    # Initialise the length of the curve as zero.
    length = 0.0

    # For each line segment joining neighbouring nodes...
    for i in xrange(number_of_inner_nodes+1):
        # Add the length of the line segment - as approximated by the trapezoidal quadrature rule to the total length
        length += 0.5 * la.mass_norm(np.subtract(curve[i + 1], curve[i]), mass_vector) * metric[i+1][0]

    # Return the total approximate length of the curve.
    return length


def GradLength(curve, metric, number_of_inner_nodes, mass_vector):
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

    # Initialise the gradient of the length as an empty list.
    grad = []

    # For each line segment joining neighbouring nodes...
    for i in xrange(1, number_of_inner_nodes+1):
        # Pre-compute the vector curve[i+1]-curve[i] to speed up the calculation.
        d_1 = np.subtract(curve[i + 1], curve[i])

        # Pre-compute the norm of curve[i+1]-curve[i] to speed up the calculation.
        d_ip1_i = la.mass_norm(d_1, mass_vector)

        # Compute a component of the gradient of the length
        t_grad = 0.5 * (d_ip1_i * (metric[i][1]/metric[i][0]) - metric[i][0] * d_1 / d_ip1_i + metric[i-1][0] *
                        np.subtract(curve[i], curve[i - 1]) / la.mass_norm(np.subtract(curve[i], curve[i - 1]),
                                                                        mass_vector))

        # Append these values into the grad list
        for t in t_grad[1:]:
            grad.append(t)

    # Return the grad variable as a NumPy array of float type.
    return np.asarray(grad, dtype='f')