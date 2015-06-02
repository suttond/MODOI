import numpy
from ase.io import read
import numpy as np


def convert_atoms_to_vector(atoms):
    """ Converts an Atomistic Simulation Environment atoms object into a vector which can be used with the curve
    shortening algorithm.

    Args:
      atoms (ase.atoms): The ASE atoms object whose position is to be converted to a vector.

    Returns:
      numpy.array: A vector of atomistic positions.

    """

    # Create a new NumPy float array.
    vector = np.asarray([], dtype='float64')

    # For each atom in the atoms object stack it's position onto the vector object.
    for atom in atoms:
        vector = np.hstack((vector, np.asarray(atom, dtype='float64')))

    # Return the resulting vector.
    return vector


def convert_vector_to_atoms(vector, dimension=3):
    """ The inverse of convert_atoms_to_vector.

    Args:
      vector (numpy.array): The vector containing the position of the atoms to be converted.
      dimension (optional int): The dimension of the space in which the simulation is taking place. The default value is 3.

    Returns:
      numpy.array: An ASE atoms object friendly NumPy array containing the atomistic positions.

    """
    return np.asarray(vector, dtype='float64').reshape((len(vector)/dimension, dimension))


def read_configuration_file(configuration_file):
    """ Import the configuration file and store the parsed values in a dictionary.

    Args:
      configuration_file (str): The address of the configuration file as required by the Python open command.

    Returns:
      dict: A dictionary containing the parsed values from the configuration file.

    """

    # Open the configuration_file in read mode.
    f = open(configuration_file, 'r')

    # For each line in the configuration file...
    for line in f.readlines():

        # Extract the 'command' which is the two letter code to the left of the equals
        command = line[:2]

        # Extract the 'value' which is the remaining text to the right of the equals, removing any whitespace
        value = line.replace(" ", "")[3:]

        # Based on the value of the command code...
        if command == 'st':
            # The location of the .xyz file describing the start configuration is given here and stored as an ASE atoms
            # object.
            molecule = read(value.strip())

            # The starting point in R^N is then extracted
            start_point = convert_atoms_to_vector(molecule.get_positions())

            # The dimension and codimension of the problem are then determined.
            dimension = len(start_point)
            codimension = dimension - 1

        elif command == 'en':
            # The location of the .xyz file describing the end configuration is given here and stored as an ASE atoms
            # object.
            end_point = convert_atoms_to_vector(read(value.strip()).get_positions())

        elif command == 'ln':
            # The parameter describing the local number of nodes used in the geodesic computation.
            local_num_nodes = int(value)

        elif command == 'gn':
            # The parameter describing the global number of nodes used in the geodesic computation.
            global_num_nodes = int(value)

        elif command == 'pa':
            # The parameter which gives the value for E as required by the Maupertuis principle.
            metric_parameters = numpy.fromstring(value, dtype='f', sep=',')

        elif command == 'to':
            # The tolerance after which we stop running the global algorithm.
            tol = float(value)

    # Return the dictionary.
    return {
        'start_point': start_point,
        'end_point': end_point,
        'dimension': dimension,
        'codimension': codimension,
        'local_number_of_nodes': local_num_nodes,
        'global_number_of_nodes': global_num_nodes,
        'metric_parameters': metric_parameters,
        'tolerance': tol,
        'molecule': molecule
    }
