from multiprocessing.connection import Client
import time
import logging
import socket

import numpy as np

from SimulationUtilities import Configuration_Processing
from SimulationUtilities.Communication_Codes import comm_code
import LinearAlgebra as la
from CustomBFGS import find_geodesic_midpoint
from MetricValues import shutdown_metric


class SimulationClient:
    """

    The purpose of this object is to compute local geodesics using a modified BFGS method. The object receives a pair of
    end points to compute the local geodesic between. The simulation client then returns a new position for the node.
    The simulation client needs the value of the potential and it's gradient function, in order to achieve this it makes
    calls to it's assigned SimulationPotential servers.

    Attributes:
      CONFIGURATION (dict) :
          A dictionary containing the parsed values from the file in configuration_file.
      CURVE_ADDRESS (str, int) :
          A tuple containing a string representing the hostname/IP and an integer for the running SimulationServer.
      AUTHKEY (str) :
          A string containing the authorisation key for the listener method.
      DELAY (float) :
          The length of time the SimulationClient should wait if there is no new available jobs, before attempting to
          contact the SimulationServer again.
      METRIC_SERVERS :
          A list containing tuples of addresses for SimulationPotential instances.
      ID (str) :
          A string that uniquely identifies the client amongst all other clients in the computation.
      MASS_MATRIX (numpy.array) :
          A NumPy matrix containing the mass matrix of the molecular system. Produced automatically from the Atomistic
          Simulation Environment.

    """
    def __init__(self, simulation_client_id, server_host, server_port, authkey, metric_server_addresses,
                 configuration_file, logfile=None, log_level=logging.INFO, callback_delay=1.0):
        """The constructor for the SimulationClient class.

        Note:
          This class is intended to be used in conjunction with running SimulationServer and SimulationPotential
          objects. It will cause a runtime exception if this condition isn't satisfied.

        Args:
          simulation_client_id (str) :
              A string that uniquely identifies the client amongst all other clients in the computation.
          server_host (str) :
              The TCP/IP hostname of the running SimulationServer instance.
          server_port (int) :
              The port number that the SimulationServer instance is communicating on.
          authkey (str, optional) :
              Authentication key used to secure process communications. Default to None for local computations to
              increase speed.
          metric_server_addresses :
              A list containing tuples of the type (str, int) containing the hostnames and ports for the running
              SimulationPotential instances.
          configuration_file (str) :
              Directory and filename of the configuration file.
          logfile (str, optional) :
              Directory and filename of the log file. Is created if doesn't exist, overwritten if it does.
          log_level (int, optional) :
              Specify level of logging required as described in the logging package documentation.
          callback_delay (float) :
              The length of time the SimulationClient should wait if there is no new available jobs, before attempting
               to contact the SimulationServer again.

        """
        # Set the SimulationClient log output to write to logfile at prescribed log level if specified. Otherwise write
        # to console output. Setting to DEBUG will cause poor performance and should only be used to debug.
        if logfile is not None:
            logging.basicConfig(filename=logfile, level=log_level, filemode='w')
        else:
            logging.basicConfig(level=logging.INFO)

        # Read configuration from configuration_file and store in SimulationPotential's CONFIGURATION attribute.
        self.CONFIGURATION = Configuration_Processing.read_configuration_file(configuration_file)

        # Set ADDRESS and AUTHKEY attributes for Client object in the start_client method.
        self.CURVE_ADDRESS = (server_host, server_port)
        self.AUTHKEY = authkey

        # Set the callback delay as described in the attributes.
        self.DELAY = callback_delay

        # Store the ADDRESS and AUTHKEY attributes for Client objects in the start_client method used to compute the
        # metric values.
        self.METRIC_SERVERS = metric_server_addresses

        # Set the client's unique identifier.
        self.ID = simulation_client_id

        # Compute the mass matrix for the molecular system.
        self.MASS_MATRIX = np.diag(np.dstack((self.CONFIGURATION['molecule'].get_masses(),) *
                                             (self.CONFIGURATION['dimension'] /
                                              len(self.CONFIGURATION['molecule'].get_masses()))).flatten())


    def start_client(self):
        """Start the instance of SimulationClient and begin computing local geodesics.

        """

        # Define a flag to indicate if contact with the SimulationServer instance is possible.
        connection_made = False

        # Create a response to send to the SimulationServer indicating that this is the first time this SimulationClient
        # has attempted to get a task.
        client_response = {'status_code': comm_code('CLIENT_FIRST_CONTACT'),
                           'client_name': self.ID}

        # Attempt to connect to the SimulationServer instance.
        try:
            # Create a Client object that communicates with the listener on CURVE_ADDRESS using password AUTHKEY.
            server = Client(self.CURVE_ADDRESS, authkey=self.AUTHKEY)

            # When a connection is made send the client message.
            server.send(client_response)

            # The client assumes the server will respond with a message, either a local geodesic to compute or a message
            # asking the client to try again after DELAY seconds.
            server_response = server.recv()

            # Interpret the servers response by first extracting the status_code variable from the response.
            server_response_code = server_response['status_code']

            # Close the connection to the server at this point to allow other clients to communicate with the
            # SimulationServer.
            server.close()

            # Store in the connection_made flag that it was possible to create a connection.
            connection_made = True

        # If it isn't possible to connect to the server than a socket.error exception is raised.
        except socket.error:
            # Write an error to the log for this client indicating that the connection couldn't be made.
            logging.warning('Failed to Make Connection to SimulationServer. Shutting down client.')

            # Send a signal to the running instances of SimulationPotential that the SimulationClient would have used
            # indicating that they should also shutdown.
            shutdown_metric(self.METRIC_SERVERS, self.AUTHKEY)

        # This is the main loop of the SimulationClient - the program stops running when it is no longer possible to
        # communicate with the SimulationServer. This is decided by the connection_made flag.
        while connection_made:

            # At this point in the code a new server_response should have been received. How the SimulationClient reacts
            # depends on the communication code received.

            # If the server has indicated it is giving the SimulationClient a new geodesic to compute then...
            if server_response_code == comm_code('SERVER_GIVES_NEW_TASK'):

                # Compute the rescaled tangent direction of the curve as store as a NumPy array.
                tangent_direction = (1 / float(self.CONFIGURATION['local_number_of_nodes'] + 1)) * \
                    np.subtract(server_response['right_end_point'], server_response['left_end_point'], dtype='float64')

                # Compute the local geodesic using the BFGS method and store the NumPy array in result
                result = \
                    find_geodesic_midpoint(server_response['left_end_point'],
                                                server_response['right_end_point'],
                                                self.CONFIGURATION['local_number_of_nodes'],
                                                la.orthonormal_tangent_basis(tangent_direction,
                                                                          self.CONFIGURATION['dimension']),
                                                tangent_direction, self.CONFIGURATION['codimension'],
                                                self.METRIC_SERVERS,
                                                self.MASS_MATRIX,
                                                self.AUTHKEY)

                # If the function find_geodesic_midpoint returned a None object then it couldn't contact it's
                # SimulationPotential instances and should be restarted.
                if result is None:
                    # Tell the user via the log that the SimulationPotential instances couldn't be contacted.
                    logging.warning('Failed to Make Connection to SimulationPotential. Shutting down client.')

                    # Exit the main loop of the SimulationClient.
                    break

                # If there is a midpoint then construct a client response to tell the server which node has which new
                # position.
                client_response = {'status_code': comm_code('CLIENT_HAS_MIDPOINT_DATA'),
                                   'node_number': server_response['node_number'],
                                   'new_node_position': result,
                                   'client_name': self.ID
                                   }

            # Otherwise if the server has asked the SimulationClient to try again later...
            elif server_response_code == comm_code('SERVER_REQUEST_CALLBACK'):
                # Make the SimulationClient wait for DELAY seconds
                time.sleep(self.DELAY)

                # Create a response to tell the SimulationServer that the SimulationClient would like a new job.
                client_response = {'status_code': comm_code('CLIENT_HAS_NO_TASK'), 'client_name': self.ID}

            # Attempt to connect to the SimulationServer instance.
            try:
                # Create a Client object that communicates with the listener on CURVE_ADDRESS using password AUTHKEY.
                server = Client(self.CURVE_ADDRESS, authkey=self.AUTHKEY)

                # When a connection is made send the client message.
                server.send(client_response)

                # The client assumes the server will respond with a message, either a local geodesic to compute or a
                # message asking the client to try again after DELAY seconds.
                server_response = server.recv()

                # Interpret the servers response by first extracting the status_code variable from the response.
                server_response_code = server_response['status_code']

                # Close the connection to the server at this point to allow other clients to communicate with the
                # SimulationServer.
                server.close()

            # If it isn't possible to connect to the server than a socket.error exception is raised.
            except (socket.error, EOFError):
                # Write an error to the log for this client indicating that the connection couldn't be made.
                logging.warning('Failed to Make Connection to SimulationServer. Shutting down client.')

                # Send a signal to the running instances of SimulationPotential that the SimulationClient would have
                # used indicating that they should also shutdown.
                shutdown_metric(self.METRIC_SERVERS, self.AUTHKEY)

                # Exit the main loop of the SimulationClient.
                break