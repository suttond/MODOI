from multiprocessing.connection import Listener
import pickle
import logging
import time

from SimulationUtilities.Configuration_Processing import read_configuration_file
from SimulationUtilities.Communication_Codes import comm_code
from SimulationUtilities.Curve import Curve
from SimulationUtilities.Visualize import write_xyz_animation


class SimulationServer:
    """

    The purpose of this object is to manage the positions of the global nodes. The object determines the next node
    to be sent using the next() method of the Curve class. The node is then sent using TCP via the Listener object
    provided by the multiprocessing package. When the curve movement drops below a tolerance specified in the
    configuration file the server will stop.

    Attributes:
      CONFIGURATION (dict) :
          A dictionary containing the parsed values from the file in configuration_file.
      ADDRESS (str, int) :
          A tuple containing a string representing the hostname/IP and an integer for the service port.
      AUTHKEY (str) :
          A string containing the authorisation key for the listener method.
      OUTPUT_FILENAME (str) :
          The directory and filename for the output. Note that this should not have a suffix as this is provided at
          runtime.
      CURVE (Curve) :
          A curve object that contains the result of the curve shortening procedure.
      FINISHED (bool) :
          To indicate whether run_simulation has completed.

    """
    def __init__(self, configuration_file, output_filename, logfile=None,
                 hostname='localhost', port=5000, authkey=None, log_level=logging.DEBUG, timeout=5000):
        """The constructor for the SimulationServer class.

        Note:
          This class is intended to be used in conjunction with running SimulationClient and SimulationPotential
          objects. See the run_simulation method for more details.

        Args:
          configuration_file (str) :
              Directory and filename of the configuration file.
          output_filename (str) :
              Directory and filename of the output files. A suffix is added at runtime and is not required.
          logfile (str, optional) :
              Directory and filename of the log file. Is created if doesn't exist, overwritten if it does.
          hostname (str, optional) :
              Hostname/IP which the SimulationServer will run on. Default value is 'localhost' to run simulations
              locally.
          port (int, optional) :
              Port number which the SimulationServer will run on. Default value is 5000, just because.
          authkey (str, optional) :
              Authentication key used to secure process communications. Default to None for local computations to
              increase speed.
          log_level (int, optional) :
              Specify level of logging required as described in the logging package documentation.
        """

        # Set the SimulationServer log output to write to logfile at prescribed log level if specified. Otherwise write
        # to console output.
        if logfile is not None:
            logging.basicConfig(filename=logfile, level=log_level, filemode='w')
        else:
            logging.basicConfig(level=logging.INFO)

        # Read configuration from configuration_file and store in SimulationServer's CONFIGURATION attribute.
        self.CONFIGURATION = read_configuration_file(configuration_file)

        # Add an additional property to the CONFIGURATION attribute containing a calculation of the total number of
        # nodes. This is used to determine the average movement of the curve.
        self.CONFIGURATION['total_number_of_nodes'] = self.CONFIGURATION['global_number_of_nodes'] * \
            (self.CONFIGURATION['local_number_of_nodes'] - 1) + 1

        # Set ADDRESS and AUTHKEY attributes for Listener object in the run_simulation method.
        self.ADDRESS = (hostname, port)
        self.AUTHKEY = authkey

        # Set OUTPUT_FILENAME attribute for save_simulation method.
        self.OUTPUT_FILENAME = output_filename

        # Create a new Curve object, default state is a straight line joining the start point to the end point
        self.CURVE = Curve(self.CONFIGURATION['start_point'], self.CONFIGURATION['end_point'],
                                 self.CONFIGURATION['global_number_of_nodes'],
                                 self.CONFIGURATION['total_number_of_nodes'])

        # Add configuration from configuration_file to the CURVE attribute. This is so that researchers can share
        # the results from the save_simulation method and determine for what parameters the simulation was run under.
        self.CURVE.configuration = self.CONFIGURATION

        # Set flag to indicate if simulation is complete
        self.FINISHED = False

        self.TIMEOUT = timeout

    def run_simulation(self):
        """Start the SimulationServer listener and start the Birkhoff curve shortening procedure.

        Note:
          This class assumes that after starting there will be at least one running instance of SimulationClient
          pointing to ADDRESS, otherwise this process will remain indefinitely blocked.
        """

        # Set up the listener for communication at ADDRESS. Backlog is set to 5 to allow maximum concurrent connections.
        # If a larger scale project is needed then the SimulationServer will need to run a separate communication
        # handling process.
        logging.info('Starting Server on %s', str(self.ADDRESS))
        server = Listener(self.ADDRESS, authkey=self.AUTHKEY, backlog=5)

        client_monitor = {}

        while True:

            # The Simulation server only receives requests from running instances of SimulationServer objects. The
            # SimulationServer waits for an incoming connection and is blocked until it receives one.
            logging.debug('Listening for connection from Client instance...')
            client = server.accept()

            # Receive request from a SimulationClient. Request is in form of dictionary that is stored in
            # client_response
            logging.debug('Connection to Client received from %s.', server.last_accepted)
            client_response = client.recv()

            # If the SimulationClient identifies that it contains a new midpoint position then extract it and update
            # the CURVE attribute.
            if client_response['status_code'] == comm_code('CLIENT_HAS_MIDPOINT_DATA'):
                logging.debug('Client response contains new midpoint.')
                # If the client hasn't already been presumed dead then update node position, otherwise ignore.
                if client_monitor[client_response['client_name']] is not None:
                    self.CURVE.set_node_position(client_response['node_number'], client_response['new_node_position'])
                else:
                    logging.debug('Client took too long to respond. Was presumed dead.')
            elif client_response['status_code'] == comm_code('CLIENT_FIRST_CONTACT'):
                logging.debug('First contact from Client:' + str(client_response['client_name']))
                client_monitor[client_response['client_name']] = None

            # Check whether each node in the global curve has now been repositioned. If so, check the total movement
            # of the curve, if this is smaller than the tolerance set in configuration_file then exit the main loop.
            # Otherwise, set all of the nodes in the global curve as movable again. For more information see the Curve
            # class.
            if self.CURVE.all_nodes_moved():
                logging.info('Total curve movement: %s', self.CURVE.movement)
                if self.CURVE.movement < self.CONFIGURATION['tolerance']:
                    break

                self.CURVE.set_node_movable()

            # If the code has reached this point then there are still nodes to move, and the desired solution hasn't
            # yet been attained. Use the next_movable_node method of the CURVE to get the index of the next movable
            # node.
            next_node_number = None
            for client_info in client_monitor:
                if client_monitor[client_info] is not None:
                    time_elapsed_since_sent = time.time() - client_monitor[client_info][0]
                    if (time_elapsed_since_sent) > self.TIMEOUT:
                        next_node_number = client_monitor[client_info][1]
                        client_monitor[client_info] = None
            if next_node_number is None:
                next_node_number = self.CURVE.next_movable_node()

            logging.debug('Next movable node: %s', next_node_number)

            # Even if there are nodes that need to be tested, it may not be possible if it's neighbours are currently
            # being tested, in which case the SimulationServer tells the client to try again later. Otherwise a new
            # node is obtained and sent back to the server.
            if next_node_number is not None:
                logging.debug('Sending new node to Client.')
                server_response = {'status_code': comm_code('SERVER_GIVES_NEW_TASK'),
                                   'node_number': next_node_number,
                                   'left_end_point': self.CURVE.get_points()[next_node_number - 1],
                                   'right_end_point': self.CURVE.get_points()[next_node_number + 1]
                                   }
                client.send(server_response)
                client_monitor[str(client_response['client_name'])] = [time.time(), next_node_number]
            else:
                logging.debug('No node available to move. Requesting callback.')
                client.send({'status_code': comm_code('SERVER_REQUEST_CALLBACK')})

            # Close connection to client
            client.close()

        # The computation has now been completed and the Listener is shut down.
        logging.info('Shutting down Server.')
        server.close()

        self.FINISHED = True

    def save_simulation(self):
        """ Save the results of the simulation to a pickle file and XYZ animation.

        """

        # Check if the simulation has finished.
        if self.FINISHED is True:
            # Save the Curve object in CURVE as a pickled object in the location specified in OUTPUT_FILENAME
            pickle.dump(self.CURVE, open(self.OUTPUT_FILENAME + '.pkl', "wb"))

            # Use the custom function write_xyz_animation to produce an XYZ file output containing an animation of the
            # trajectory.
            write_xyz_animation(self.OUTPUT_FILENAME + '.pkl', self.OUTPUT_FILENAME + '.xyz')
        else:
            logging.warn('save_simulation failed as simulation has not been run.')

    def set_output_file(self, output_filename):
        """ Mutator method for ad hoc change of output filename.

        """
        if self.FINISHED is True:
            self.OUTPUT_FILENAME = output_filename
        else:
            logging.warn('Output filename not changed as method accessed before simulation finished.')

    def get_curve(self):
        """ Accessor method for  Curve.

        """
        if self.FINISHED is True:
            return self.CURVE
        else:
            logging.warn('No curve object returned as simulation is not finished.')
            return None