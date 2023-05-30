from classes import Molecule

# Parses a file with one or more molecules in mol2 format
def parse_mol2_file(file_path):
    
    """
    Parses a file with one or more molecules in mol2 format.

    Args:
        file_path (str): The path to the .mol2 file to parse.

    Returns:
        list: A list of Molecule objects parsed from the .mol2 file.
    """

    with open(file_path, 'r') as file:
        content = file.read()
    molecules_str = content.split('@<TRIPOS>MOLECULE')[1:]
    molecules = [Molecule('@<TRIPOS>MOLECULE' + mol_str) for mol_str in molecules_str]
    return molecules