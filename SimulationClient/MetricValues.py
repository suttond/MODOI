import math
import socket
import logging
from multiprocessing.connection import Client

from SimulationUtilities.Communication_Codes import comm_code


def get_metric(curve, number_of_inner_points, metric_server_addresses, authkey):
    """ This function distributes the task of computing the metric values along the curve using the SimulationPotential instances.

    Args:
      curve (numpy.array): A list of NumPy arrays representing a local geodesic.
      number_of_inner_points (int): The number of points along the curve, less two.
      metric_server_addresses: A list of tuples of the form (str, int) containing the hostnames and port numbers of the SimulationPotential instances.
      authkey (str): The password used in order to communicate with the SimulationPotential instances.

    Returns:
      list: A list of float values called metric where metric[i] = a(curve[i]).

    """

    # Initialise the memory for the metric values
    metric = [[]] * (number_of_inner_points + 2)

    # Compute how many SimulationPotential instances are available to the SimulationClient
    number_of_metric_servers = len(metric_server_addresses)

    # Estimate how many tasks per SimulationPotential is a uniform distribution
    tasks_per_server = int(math.floor((number_of_inner_points + 2) / number_of_metric_servers))

    # Compute how many tasks are left over if I give each SimulationPotential tasks_per_server points to evaluate.
    remainder_tasks = number_of_inner_points + 2 - tasks_per_server * number_of_metric_servers

    # Create a flag to indicate which point along the curve is currently being issued
    current_point = 0

    # Indicate whether the current PotentialServer instance needs an extra task to ensure all the nodes are computed
    extra_jobs = 1

    # For each SimulationPotential instance...
    for address in xrange(number_of_metric_servers):

        # If all of the surplus tasks have been issued then indicate no extra tasks are needed
        if address >= remainder_tasks:
            extra_jobs = 0

        # Initialise a list to contain the points to be computed for the current SimulationPotential instance
        points = []

        # Get the points between and including curve[current_point] to curve[current_point+tasks_per_server+extra_jobs]
        # and prepare them to send to the current SimulationPotential instance
        for i in xrange(current_point, current_point + tasks_per_server + extra_jobs):
            points.append([curve[i], i])

        # Increase current_point so the same points are not sent again
        current_point += tasks_per_server + extra_jobs

        # Attempt to contact the SimulationPotential instance
        try:
            # Create connection with the SimulationPotential instance
            metric_server = Client(metric_server_addresses[address], authkey=authkey)

            # Prepare a response containing the points to evaluate the metric on
            metric_server.send({'status_code': comm_code('CLIENT_PROVIDES_POINT'),
                                'points': points})

            # Close the connection
            metric_server.close()

        # If the SimulationPotential instance couldn't be contacted then return None to shut down the SimulationClient
        except socket.error:
            # Write a warning to the log explaining which SimulationPotential couldn't be contacted.
            logging.warning('Failed to Make Connection to SimulationPotential at '
                          + str(metric_server_addresses[address]) + '.')

            # Return a None value
            return None

    # Now the code attempts to collect the computed values.

    # For each SimulationPotential instance...
    for address in xrange(number_of_metric_servers):
        # Attempt to contact the SimulationPotential instance
        try:
            # Create connection with the SimulationPotential instance
            metric_server = Client(metric_server_addresses[address], authkey=authkey)

            # Prepare a response indicating the SimulationClient is ready for the answer
            metric_server.send({'status_code': comm_code('CLIENT_ASKS_FOR_VALUES')})

            # Receive the metric values at the given points
            metric_server_response = metric_server.recv()

            # Process the received values into the metric list
            for value in metric_server_response['values']:
                metric[value[1]] = value[0]

            # Close the connection
            metric_server.close()

        # If the SimulationPotential instance couldn't be contacted then return None to shut down the SimulationClient
        except socket.error:
            # Write a warning to the log explaining which SimulationPotential couldn't be contacted.
            logging.warning('Failed to Make Connection to SimulationPotential at '
                          + str(metric_server_addresses[address]) + '.')

            # Return a None value
            return None

    # If None hasn't been returned then return the metric values
    return metric


def shutdown_metric(metric_server_addresses, authkey):
    """ This function tells all of the SimulationPotential instances running on the addresses in metric_server_addresses to shutdown. This is called when a SimulationClient instance shuts down.

    Args:
      metric_server_addresses: A list of tuples of the form (str, int) containing the hostnames and port numbers of the SimulationPotential instances.
      authkey (str): The password used in order to communicate with the SimulationPotential instances.

    Returns:
      float: The length of the curve with metric values in metric.

    """

    # For each SimulationPotential instance...
    for address in xrange(len(metric_server_addresses)):

        # Try making contact with the SimulationPotential instance...
        try:
            # Establish a connection with the SimulationPotential
            metric_server = Client(metric_server_addresses[address], authkey=authkey)

            # Send a status code indicating the SimulationPotential instance should stop running.
            metric_server.send({'status_code': comm_code('KILL')})

            # Close the connection.
            metric_server.close()

        # If contact with the SimulationPotential instance cannot be made then...
        except socket.error:
            # Make a note in the log which SimulationPotential couldn't be contacted.
            logging.warning('Failed to Make Connection to SimulationPotential at '
                          + str(metric_server_addresses[address]) + '.')