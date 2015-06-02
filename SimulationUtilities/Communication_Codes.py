

def comm_code(label):
    """ This function takes a string describing the type of communication being made and converts it to a machine
    friendly value.

    Args:
      label (str): The human-readable form of the communication code.

    Returns:
      int: The machine-readable form of the communication code.

    """

    # Create a dictionary mapping the communication codes to numbers
    codes = {
            'CLIENT_HAS_NO_TASK': 0,
            'CLIENT_HAS_MIDPOINT_DATA': 1,
            'SERVER_GIVES_NEW_TASK': 2,
            'SERVER_REQUEST_CALLBACK': 3,
            'CLIENT_PROVIDES_POINT': 4,
            'SERVER_PROVIDES_VALUES': 5,
            'CLIENT_ASKS_FOR_VALUES': 6,
            'CLIENT_HEARTBEAT': 7,
            'KILL': 8,
            'CLIENT_FIRST_CONTACT': 9
    }

    # Try sending back the numeric value for the label inputted
    try:
        return codes[label]
    # If the wrong label is requested return a None object.
    except KeyError:
        return None