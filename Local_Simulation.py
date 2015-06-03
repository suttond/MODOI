import multiprocessing
import logging

from SimulationServer.SimulationServer import SimulationServer
from SimulationClient.SimulationClient import SimulationClient
from SimulationPotential.SimulationPotential import SimulationPotential


def startMdServer(configuration_file, output_filename, logfile, hostname, port, authkey, loglevel):
    """ A function to create and start a SimulationServer instance.

    """
    md_server = SimulationServer(configuration_file, output_filename, logfile, hostname, port, authkey, loglevel)
    md_server.run_simulation()
    md_server.save_simulation()


def startMetServer(configuration_file, log_file, log_level, hostname, port, authkey):
    """ A function to create and start a SimulationPotential instance.

    """
    met_server = SimulationPotential(configuration_file, log_file, log_level, hostname, port, authkey)
    met_server.run_potential_server()


def startMdClient(client_id, server_host, server_port, authkey, metric_server_addresses, configuration_file, logfile):
    """ A function to create and start a SimulationServer instance.

    """
    md_client = SimulationClient(client_id, server_host, server_port, authkey, metric_server_addresses,
                                 configuration_file, logfile)
    md_client.start_client()

if __name__ == '__main__':

    # Configure distributed system for head machine
    NUMBER_OF_METRIC_SERVERS_PER_CLIENT = 2
    NUMBER_OF_CLIENTS = 2
    IP_ADDRESS = 'localhost'
    SERVER_PORT = 5000
    BASE_PORT = 5000  # Typically same value as server port
    AUTHKEY = 'password'

    # Location of the configuration file
    CONFIG_FILE = 'Experiment/Example.bkhf'

    # Location and name of the output file without a suffix (will be added in runtime)
    OUTPUT = 'Experiment/Trajectory'

    # Create a local process running a SimulationServer instance
    s = multiprocessing.Process(target=startMdServer, args=(CONFIG_FILE, OUTPUT, None, IP_ADDRESS, BASE_PORT,
                                                            AUTHKEY, logging.INFO))
    s.start()

    for i in xrange(NUMBER_OF_CLIENTS):
        metric_server_addresses = []

        # For the current client create NUMBER_OF_METRIC_SERVERS_PER_CLIENT many SimulationPotential instances
        for j in xrange(NUMBER_OF_METRIC_SERVERS_PER_CLIENT):
            # Compute a new port number (to avoid local conflicts)
            BASE_PORT += 1
            metric_server_addresses.append(('localhost', BASE_PORT))
            m = multiprocessing.Process(target=startMetServer, args=(CONFIG_FILE, None,
                                                                     logging.INFO, IP_ADDRESS, BASE_PORT, AUTHKEY, ))
            m.start()

        # Start a SimulationClient process in a separate thread.
        UNIQUE_CLIENT_ID = 'Client_' + str(i)
        c = multiprocessing.Process(target=startMdClient, args=(UNIQUE_CLIENT_ID, IP_ADDRESS, SERVER_PORT, AUTHKEY,
                                                                metric_server_addresses, CONFIG_FILE, None))
        c.start()