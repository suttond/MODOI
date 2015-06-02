import ase

def extract_states_pdb(pdb_filename, x0_location, xN_location, x0_index, xN_index):
    """ This function is an ASE wrapper that produces a start and end configuration from a PDB file containing multiple
    conformation.

    Args:
      pdb_filename (str): The location of the PDB file.
      x0_location (str): The location of where to write the x0.xyz file, that is the starting configuration.
      xN_location (str): The location of where to write the xN.xyz file, the final configuration.
      x0_index (int): The index in the PDB file of the starting configuration.
      xN_index (int): The index in the PDB file of the final configuration.

    """
    ase.io.write(x0_location, ase.io.read(pdb_filename, index=x0_index))
    ase.io.write(xN_location, ase.io.read(pdb_filename, index=xN_index))