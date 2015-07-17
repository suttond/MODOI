from multiprocessing.connection import Listener
import math
import logging

from ase.calculators.emt import EMT

from SimulationUtilities import Configuration_Processing
from SimulationUtilities.Communication_Codes import comm_code


class SimulationPotential:
    """

    The purpose of this object is to provide an independent service that computes the forces and potential energy for
    the system. The object listens for a TCP connection containing a molecular configuration. The forces and potential
    energy are then sent using TCP via the Listener object provided by the multiprocessing package. If a kill code is
    received then the server shuts down.

    Attributes:
      CONFIGURATION (dict) :
          A dictionary containing the parsed values from the file in configuration_file.
      ADDRESS (str, int) :
          A tuple containing a string representing the hostname/IP and an integer for the service port.
      AUTHKEY (str) :
          A string containing the authorisation key for the listener method.


    """
    def __init__(self, configuration_file, logfile=None, log_level=logging.INFO,
                 hostname='localhost', port=5001, authkey='password'):
        """The constructor for the SimulationPotential class.

        Note:
          This class is intended to be used in conjunction with running SimulationServer and SimulationClient objects.
          See the run_potential_server method for more details.

        Args:
          configuration_file (str) :
              Directory and filename of the configuration file.
          logfile (str, optional) :
              Directory and filename of the log file. Is created if doesn't exist, overwritten if it does.
          log_level (int, optional) :
              Specify level of logging required as described in the logging package documentation.
          hostname (str, optional) :
              Hostname/IP which the SimulationPotential will run on. Default value is 'localhost' to run simulations
              locally.
          port (int, optional) :
              Port number which the SimulationPotential will run on. Default value is 5000, just because.
          authkey (str, optional) :
              Authentication key used to secure process communications. Default to None for local computations to
              increase speed.

        """

        # Set the SimulationPotential log output to write to logfile at prescribed log level if specified. Otherwise
        # write to console output. Setting to DEBUG will cause poor performance and should only be used to debug.
        if logfile is not None:
            logging.basicConfig(filename=logfile, level=log_level, filemode='w')
        else:
            logging.basicConfig(level=logging.INFO)

        # Read configuration from configuration_file and store in SimulationPotential's CONFIGURATION attribute.
        self.CONFIGURATION = Configuration_Processing.read_configuration_file(configuration_file)

        # Set ADDRESS and AUTHKEY attributes for Listener object in the run_simulation method.
        self.ADDRESS = (hostname, port)
        self.AUTHKEY = authkey

    def run_potential_server(self, small_number=1e-12):
        """Start the instance of SimulationPotential ready to receive requests for metric values.

        Note:
          This class assumes that after starting there will be at least one running instance of SimulationClient
          pointing to ADDRESS, otherwise this process will remain indefinitely blocked.

        Args:
          small_number (str) :
              A small number used to represent the zero metric value. This is used instead of zero to prevent divide
              by zero errors when evaluating the gradient of the metric.

        """

        # Extract the ASE atoms object for molecule and set the calculator to the pure Python EMT implementation
        molecule = self.CONFIGURATION['molecule']
        molecule.set_calculator(EMT())

        # Set up the listener for communication at ADDRESS. This doesn't queue connections as it is **assumed** that
        # each SimulationClient has it's own pool of SimulationPotential servers.
        logging.info('Starting Potential Server on %s', str(self.ADDRESS))
        server = Listener(self.ADDRESS, authkey=self.AUTHKEY)

        while True:

            # The Simulation server only receives requests from running instances of SimulationServer objects. The
            # SimulationServer waits for an incoming connection and is blocked until it receives one.
            logging.debug('Listening for connection from SimulationClient instance...')
            client = server.accept()

            # Receive request from a SimulationClient. Request is in form of dictionary that is stored in
            # client_response
            logging.debug('Connection to SimulationClient received from %s.', server.last_accepted)
            client_response = client.recv()

            # Prepare/reset list container for server response
            values = []

            # This is the main decision logic for handling a connection from a client server.
            # If the SimulationClient indicates it is providing points to evaluate the metric on then compute those
            # values and prepare a response.
            if client_response['status_code'] == comm_code('CLIENT_PROVIDES_POINT'):
                logging.debug('Client provides point data to evaluate.')

                # Close connection before computing metric values. To free up network traffic.
                client.close()

                # For each point received in client request, update the positions in our ASE atoms object and then
                # compute the potential energy and forces.
                for point in client_response['points']:
                    molecule.set_positions(Configuration_Processing.convert_vector_to_atoms(point[0]))
                    values.append([[math.sqrt(max([self.CONFIGURATION['metric_parameters'][0] -
                                                   molecule.get_potential_energy(), small_number])),
                                   Configuration_Processing.convert_atoms_to_vector(molecule.get_forces())], point[1]])

                # Compile response ready for when client asks for it.
                server_response = {'status_code': comm_code('SERVER_PROVIDES_VALUES'), 'values': values}

            elif client_response['status_code'] == comm_code('CLIENT_ASKS_FOR_VALUES'):

                # Send computed potential energies and forces to client.
                logging.debug('Client requests result from computation.')
                client.send(server_response)
                client.close()

            elif client_response['status_code'] == comm_code('KILL'):

                # Break main loop of SimulationPotential and subsequently close the server.
                logging.debug('Shut-down signal received.')
                break

        # Close the SimulationPotential
        logging.info('Shutting down SimulationPotential.')
        server.close()