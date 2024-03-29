#This file contains the definition for the Molecule class which is used to represent and process molecular data

import pandas as pd
from pymongo import MongoClient, errors
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, rdmolops
from openbabel import openbabel
from openbabel import pybel
import pkasolver
from pkasolver.query import calculate_microstate_pka_values
import logging
import pickle
from datetime import datetime
import math
import requests
import os
import json
import time
import glob

import os
import tempfile
import subprocess
from subprocess import Popen, PIPE, STDOUT
from shutil import copyfile

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

def compute_scaffold(smiles):
    """
    Compute the Bemis-Murcko scaffold of a molecule given its SMILES representation.

    Parameters
    ----------
    smiles : str
        The SMILES string of the molecule.

    Returns
    -------
    scaffold_smiles : str
        The SMILES string of the computed scaffold. If an error occurs during computation, 
        this function returns math.nan.

    Example
    -------
    compute_scaffold(smiles="CCO")

    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smiles = Chem.MolToSmiles(scaffold)
        return scaffold_smiles
    except Exception as e:
        print(f"Failed to compute scaffold for SMILES {smiles}: {e}")
        return math.nan

class Molecule:
    """
    The Molecule class is designed to represent and process molecular data. On a larger scale, the 'Molecule'
    class handles single-item operations, while the 'Database' class handles many-item operations. They both
    work together to manage the mongo operations.

    Parameters
    ----------
    smiles : str
        The SMILES string representation of the molecule.
    id : str, optional
        The ID of the molecule. Defaults to None.
    input_format : str, optional
        The format of the input_data. Defaults to None.
    input_data : str, optional
        The data of the molecule in the specified input_format. Defaults to None.

    Attributes
    ----------
    mol2 : str
        The molecule in mol2 format.
    sections : dict
        A dictionary containing parsed sections from the mol2 data.

    Example
    -------
    mol = Molecule(smiles="CCO")
    """
    def __init__(self, smiles, id=None, input_format=None, input_data=None):
        self.smiles = smiles
        self.mol_id = id
        self.input_format = input_format
        self.input_data = input_data
        if self.input_data is None:
            self.mol2 = self._get_mol2()
            
        else:
            self.mol2 = self.convert_format('mol2', input_data=self.input_data)
        self.sections = self._parse_mol2()
        self.dict = self.to_dict()

    def _get_mol2(self, pHs=None):
        """
        Intermediate private method. Will be deprecated soon.
        """
        return self.convert_format('mol2', pHs)

    def _parse_mol2(self):
        """
        Private method to parse the mol2 format of the molecule.

        Returns
        -------
        dict
            A dictionary containing the parsed sections from the mol2 data.
        """
        sections = self.mol2.split('@<TRIPOS>')
        # Parse MOLECULE section
        molecule_lines = sections[1].split('\n')
        self.zinc_id = molecule_lines[1].strip()
        self.num_atoms, self.num_bonds, self.num_substructures, self.num_features, self.num_sets = map(int, molecule_lines[2].split())
        self.molecule_type = molecule_lines[3].strip()
        self.charge_type = molecule_lines[4].strip()
        if len(molecule_lines) > 5:
            self.comments = molecule_lines[5:]

        # Parse the rest of the sections
        section_dict = {}
        for section in sections[2:]:
            section_lines = section.split('\n')
            section_name = section_lines[0]
            section_data = section_lines[1:]
            if section_name in ['ATOM', 'BOND', 'SUBSTRUCTURE']:
                df = pd.DataFrame([s.split() for s in section_data if s])
                section_dict[section_name.lower()] = df.to_dict()
            else:
                section_dict[section_name.lower()] = section_data

        return section_dict

    def convert_format(self, output_format, pHs=None, input_data=None):
        """
        Convert a single molecule to a specified format. Can use RDKit to correct for different pH
        conditions.

        Parameters
        ----------
        output_format : str
            The desired output format.
        pHs : list, optional
            A list of pH values for which to generate the output. Defaults to None.
        input_data : str, optional
            The input data in pdbqt format. Defaults to None.

        Returns
        -------
        str or dict
            The molecule in the specified format. If pHs is specified, a dictionary is returned 
            with each pH as the key and the corresponding molecule in the specified format as the value.
        """
        if input_data is None:
            mol = Chem.MolFromSmiles(self.smiles)
            mol = Chem.AddHs(mol)  # Add Hs
            AllChem.EmbedMolecule(mol)  # Generate 3D coordinates
            AllChem.UFFOptimizeMolecule(mol)  # universal force field optimization
            molfile = Chem.MolToMolBlock(mol)
            pybel_mol = pybel.readstring("mol", molfile)
        else:
            pybel_mol = pybel.readstring("pdbqt", input_data)
        
        if pHs is None:
            return pybel_mol.write(output_format)
        else:
            # For pH adjusted molecules
            pH_corrected_mol2s = {}
            for pH in pHs:
                obConversion = openbabel.OBConversion()
                obConversion.SetInAndOutFormats("mol", output_format)
                obConversion.AddOption("p", openbabel.OBConversion.OUTOPTIONS, str(pH))
                pH_corrected_mol2s[pH] = obConversion.WriteString(pybel_mol.OBMol)
                obConversion.RemoveOption("p", openbabel.OBConversion.OUTOPTIONS)
                
            return pH_corrected_mol2s

    def to_dict(self):
        """
        Convert the Molecule object to a dictionary.

        Returns
        -------
        dict
            A dictionary containing the attributes of the Molecule object.
        """
        return {
            'zinc_id': self.zinc_id,
            'num_atoms': self.num_atoms,
            'num_bonds': self.num_bonds,
            'num_substructures': self.num_substructures,
            'num_features': self.num_features,
            'num_sets': self.num_sets,
            'molecule_type': self.molecule_type,
            'charge_type': self.charge_type,
            'comments': getattr(self, 'comments', None),
            'sections': self.sections
        }

class Database:
    def __init__(self, name, host = 'mongo_db', port = 27017,
                 username = 'unam', password = '12345'):
        self.name = name
        self.client = MongoClient(
                host = host,
            port = port,
            serverSelectionTimeoutMS = 3000, # 3 second timeout
            username = username,
            password = password,
        )
        self.client.server_info()  # Trigger ServerSelectionTimeoutError if cannot connect to MongoDB server
        self.DB = self.client[self.name]

    def insert_old(self, collection_name, documents):
        """
        Inserts the provided documents into the specified collection.

        Args:
        collection_name (str): The name of the collection to insert documents into.
        documents (list of dict): The documents to be inserted.
        """
        collection = self.client[self.name][collection_name]
        collection.create_index("smiles", unique=True)
        try:
            collection.insert_many(documents, ordered=False)
        except errors.BulkWriteError as e:
            #print(e.details['writeErrors'])
            ''
    def insert(self, collection_name, documents):
        """
        Inserts the provided documents into the specified collection, overwriting existing documents with the same index.

        Args:
        collection_name (str): The name of the collection to insert documents into.
        documents (list of dict): The documents to be inserted.
        """
        collection = self.client[self.name][collection_name]
        collection.create_index("index", unique=True)
        try:
            for doc in documents:
                collection.update_one({"index": doc["index"]}, {"$set": doc}, upsert=True)
        except errors.BulkWriteError as e:
            print(e.details['writeErrors'])

    def exists(self, collection_name, index):
        """
        Checks if a document with the specified index exists in the specified collection.

        Args:
        collection_name (str): The name of the collection.
        index (str): The index to check.

        Returns:
        bool: True if the document exists, False otherwise.
        """
        collection = self.client[self.name][collection_name]
        return collection.count_documents({"index": index}) > 0

    def explore(self, show=True):
        """
        Retrieves and displays information about the collections in the database.

        Args:
        show (bool, optional): Whether to print the information. Defaults to True.

        Returns:
        dict: A dictionary with collection names as keys and collection info as values.
        """
        collections_info = {}
        collections = self.DB.list_collection_names()
        if show: print("\n--- Collections Information ---")
        for collection_name in collections:
            collection = self.DB[collection_name]
            collection_info = {}
            columns = collection.find_one()
            if columns:
                for column, value in columns.items():
                    column_type = type(value).__name__
                    collection_info[column] = column_type

            index_info = collection.index_information()
            indexed_columns = [index['key'][0][0] for index in index_info.values() if index['key'][0][0] != '_id']
            collection_info['indexed_columns'] = indexed_columns if indexed_columns else 'No indexed columns'

            collections_info[collection_name] = collection_info
            if show:
                for collection, info in collections_info.items():
                    print(f"\n\nCollection: {collection}\n")
                    for column, data_type in info.items():
                        if column == 'indexed_columns':
                            print(f"{column}: {data_type}")
                        else:
                            print(f"{column}\tData Type: {data_type}")

        return collections_info

    def smiles_from_sdf(self, collection_name):
        import tqdm
        collection = self.DB[collection_name]
        temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False)
        print(f'Creating new smiles column from the sdf column in {collection_name}...')
        for document in collection.find({}):
            sdf_data = document['sdf']
            
            # Write the sdf_data to a temporary file
            temp_file.seek(0)
            temp_file.write(sdf_data)
            temp_file.truncate()
            
            # Read the SDF data from the file using RDKit
            suppl = Chem.SDMolSupplier(temp_file.name)
            mol = next(suppl)

            # sometimes rdkit fails to convert sdf files, this is a workaround
            if mol is None:
                smiles = np.nan
            else:
                smiles = Chem.MolToSmiles(mol)

            collection.update_one({"_id": document["_id"]}, {"$set": {"smiles": smiles}})
            
        temp_file.close()
        os.unlink(temp_file.name)

    def query(self, collection_name, query, store_results=True):
        """
        Executes a query on the specified collection and optionally stores the results in a new collection.

        Args:
        collection_name (str): The name of the collection to query.
        query (dict): The MongoDB query to execute.
        store_results (bool, optional): Whether to store the query results in a new collection. Defaults to True.

        Returns:
        str or Cursor: The name of the new collection (if store_results is True), or the Cursor object with the query results (if store_results is False).
        """
        collection = self.DB[collection_name]
        results = collection.find(query)
        
        if store_results:
            new_collection_name = collection_name + '_' + datetime.now().strftime('%Y%m%d%H%M%S')
            new_collection = self.DB[new_collection_name]
            new_collection.insert_many(results)
            return new_collection_name
        else:
            return results

    def convert_collection_column(self, collection_name, column_name, type_str, onError = None):
        """
        Convert an entire column in a specified collection to a certain data type.

        Parameters
        ----------
        collection_name : str
            The name of the collection in the database.
        column_name : str
            The name of the column in the collection to convert.
        type_str : str
            The type to convert the column to. For floats, pass 'double' instead of 'float'.
        onError : str, optional
            If conversion fails, original entry will be kept by default.
            To convert a float stored as a str to int, it's necessary to convert to 'double' first then to 'int'.
            By default None which keeps the original value if conversion fails.
        """
        if not onError: onError = f"${column_name}"
        collection = self.DB[collection_name]

        pipeline = [
            {
                "$addFields": {
                    column_name: {
                        "$ifNull": [
                            {
                                "$convert": {
                                    "input": f"${column_name}",
                                    "to": type_str,
                                    "onError": onError  # keep original value if conversion fails
                                }
                            },
                            f"${column_name}"  # if field does not exist, keep it as null
                        ]
                    }
                }
            },
            {"$out": collection_name},
        ]
        collection.aggregate(pipeline)

    def add_scaffold_column(self, collection, index = True):
        """
        Compute and add a scaffold column to a specified collection.

        Parameters
        ----------
        collection : str
            The name of the collection in the database.
        index : bool, optional
            Whether to create an index on the new scaffold column, by default True.
        """
        collection = self.DB[collection]
        documents = collection.find()
        
        for document in documents:
            smiles = document["smiles"]
            scaffold_smiles = compute_scaffold(smiles)
            collection.update_one({"_id": document["_id"]}, {"$set": {"scaffold": scaffold_smiles}})
        
        if index: collection.create_index("scaffold")

    def scaffold_search(self, collection_name, smiles, store_results=True):
        """
        Perform a scaffold search based on a provided SMILES string. Can create a new collection
        or act as a generator if data is large and non-persistent.

        Parameters
        ----------
        collection_name : str
            The name of the collection in the database.
        smiles : str
            The SMILES string to compute the scaffold from and search for.
        store_results : bool, optional
            Whether to store the results in a new collection, by default True.

        Returns
        -------
        str
            The name of the new collection if store_results is True.
        generator
            A generator yielding documents if store_results is False.
        """
        
        # should replace this here
        # collection.index_information()
        if 'scaffold' not in self.explore(show=False)[collection_name]:
            self.add_scaffold_column(collection_name)
        
        query_scaffold = compute_scaffold(smiles)
        if query_scaffold is None:
            print(f"Failed to compute scaffold for SMILES: {smiles}")
            return

        query = {"scaffold": query_scaffold}
        collection = self.DB[collection_name]
        results = collection.find(query)
        
        if store_results:
            new_collection_name = collection_name + '_' + datetime.now().strftime('%Y%m%d%H%M%S')
            new_collection = self.DB[new_collection_name]
            new_collection.insert_many(results)
            return new_collection_name
        else: # this else condition turns the function into a generator
            for document in results:
                yield document

    def convert_smiles_to_mol2(self, collection_name, pHs=[4, 7, 10]):
        """
        Convert SMILES strings to mol2 format for all documents in a specified collection.

        Parameters
        ----------
        collection_name : str
            The name of the collection in the database.
        pHs : list, optional
            List of pH values to use for the conversion, by default [4, 7, 10].
        """
        collection = self.DB[collection_name]
        documents = collection.find({}, {"smiles": 1})

        for document in documents:
            mol = Molecule(document["smiles"])
            mol2 = mol._get_mol2(pHs=pHs)
            
            for pH, mol2_data in mol2.items():
                column_name = f"pH_{pH}"
                
                # If the field already exists in the document, skip
                if collection.count_documents({"_id": document["_id"], column_name: {"$exists": True}}, limit = 1) > 0:
                    continue
                
                # If not, then update the document
                collection.update_one({"_id": document["_id"]}, {"$set": {column_name: mol2_data}})


    def new_format(self, collection_name, column, output_format):
        """
        Converts molecular format in a given column of a specified MongoDB collection.

        Parameters:
            collection_name: str, name of the MongoDB collection
            column: str, name of the column in the collection where the conversion is to be done
            output_format: str, the desired output molecular format (ex. 'mol2', 'sdf')
        
        Raises:
            ValueError: If the specified collection or column does not exist
        """
        # Check if collection exists
        if collection_name not in self.DB.list_collection_names():
            raise ValueError(f"Collection '{collection_name}' does not exist in the database.")

        # Check if 'smiles' field exists in the collection
        if 'smiles' not in self.explore(show=False)[collection_name]:
            raise ValueError(f"Collection '{collection_name}' does not contain a 'smiles' field.")

        collection = self.DB[collection_name]
        documents = collection.find()

        for document in documents:
            smiles = document[column]
            mol = Molecule(smiles)
            new_format_string = mol.convert_format(output_format)

            # Replace the old value in column with the new format string
            collection.update_one({"_id": document["_id"]}, {"$set": {output_format: new_format_string}})
    
    # double check logic here (!!)
    def generate_pH_corrected_mols(self, collection_name, ph_list, prepare=True):
        from meeko import MoleculePreparation
        collection = self.DB[collection_name]

        for document in collection.find({'smiles': {'$ne': float('nan')}}): # queried cursor to ignore nan values in smiles column
            smiles = document['smiles']
            print(smiles)
            print('#'*30)
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # generate a list of pkasolver state objects (protonation_states)
                protonation_states = calculate_microstate_pka_values(mol, only_dimorphite=False)
                num_protonation_states = sum(1 for _ in protonation_states)  # Count the protonation states

                # If there are no protonation states, just continue with the original molecule without pH correction.
                if num_protonation_states == 0:
                    ph_corrected_mol = mol
                else: # naive implementation of 'most common prot state at pH=x'
                    for pH in ph_list: # compare each of the molecule's pka with the working pH 
                        max_pka = -float('inf')
                        max_index = None
                        for i, state in enumerate(protonation_states):
                            if state.pka >= pH:
                                max_index = i
                                break
                            elif state.pka > max_pka:
                                max_pka = state.pka
                                max_index = i
                        # this conditional could be useful if we wanna straight up ignore molecules with no ionizable groups (would need to remove dummy counter)
                        #if max_index is None:
                            #print(f"No suitable index found for pH {pH}, SMILES: {smiles}")
                            #input('??')
                            #continue 
                    ph_corrected_mol = Chem.Mol(protonation_states[max_index].protonated_mol)
                ph_corrected_mol = Chem.AddHs(ph_corrected_mol)
                                
                # Embed conformers and store the pickled rdkit obj in mongo
                conf_ids = AllChem.EmbedMultipleConfs(ph_corrected_mol, numConfs=200, numThreads=0, pruneRmsThresh=1.5)
                pickled_mol = pickle.dumps(mol)
                phstr = str(pH*100)

                update_content = {f'pickled_rdkit_obj_pH_{phstr}': pickled_mol}

                if prepare:
                    pdbqt_strings = []
                    for conf_id in conf_ids:
                        # Make a copy of the molecule
                        mol_conf = Chem.Mol(ph_corrected_mol)

                        # Remove all conformers then read the specific one we want
                        for id in conf_ids:
                            mol_conf.RemoveConformer(id)
                        mol_conf.AddConformer(ph_corrected_mol.GetConformer(conf_id), assignId=True)

                        # Prepare the molecule with Meeko
                        preparator = MoleculePreparation()
                        preparator.prepare(mol_conf)
                        pdbqt_string = preparator.write_pdbqt_string()

                        pdbqt_strings.append(pdbqt_string)
                    print()
                    update_content[f'pdbqt_conformers_pH_{phstr}'] = pdbqt_strings
                print(f'Update data:\n\n{update_content}')
                collection.update_one({'_id': document['_id']}, {'$set': update_content})
                updated_doc = collection.find_one({'_id': document['_id']})
                if any(isinstance(val, float) and math.isnan(val) for val in updated_doc.values()):
                    print(f"NaN found after insertion for SMILES: {smiles}")

    def basic_qvina_analysis(self, collection_name):
        """
        Performs Qvina analysis on molecules in a MongoDB collection.

        Parameters:
            collection_name: str, name of the MongoDB collection

        Raises:
            ValueError: If the specified collection does not exist
        """

        '''
        
        # NOTE1: USES LOTS OF I/O
        # NOTE2: NEEDS ERROR HANDLING
 
        Takes a collection and converts its 'mol2' column into qvina outputs

  
        Assumes the contents of 'mol2' are lists of mol2 strings
        '''
        from tqdm import tqdm
        if collection_name not in self.DB.list_collection_names():
            raise ValueError(f"Collection '{collection_name}' does not exist in the database.")
    
        columns = self.explore(show=False)[collection_name]
        pH_columns = [col for col in columns if col.startswith('pdbqt_conformers_pH_')]
        collection = self.DB[collection_name]
        total_docs = collection.count_documents({})
        documents = tqdm(collection.find(), total=total_docs, desc="Processing documents")

        #collection = self.DB[collection_name]
        #documents = collection.find()
        with tempfile.NamedTemporaryFile(dir='.', suffix='.pdbqt', delete=False) as temp:
            temp_path = os.path.basename(temp.name)

        ligand_out_path = temp_path.replace('.pdbqt', '_out.pdbqt')
        tqdm.write(f'Applying qvina analysis to "{collection_name}"...')
        for doc in documents:
            tqdm.write('1')
            qvina_output_dict = {}
            for pH_col in pH_columns:
                tqdm.write('2')
                mol2_data = doc[pH_col]
                conformer_list = list(mol2_data.values())[0]
                for id, conformer in enumerate(conformer_list):
                    # Create / overwrite temp file
                    tqdm.write('3')
                    with open(temp_path, 'w') as temp:
                        temp.write(conformer)

                    #mgltools_env = "/usr/local/MGLToolsPckgs/AutoDockTools/Utilities24/"
                    orrery_env = "/opt/conda/envs/Orrery/bin/qvina2"
                    config_file = "/app/Orrery/config.vd"
                    #receptor_path = "receptor_65.pdbqt"
                    ligand_path = temp_path
                    tqdm.write('4')
                    # Naive process pipeline
                    process = subprocess.Popen(["/bin/bash"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                    command = f'{orrery_env} --config {config_file} --ligand {ligand_path}\n'
                    tqdm.write(command)
                    process.stdin.write(command.encode())

                    process.stdin.flush()
                    stdout, stderr = process.communicate()
                    tqdm.write('5')
                    ligand_out_path = temp_path.replace('.pqdbt', '_out.pdbqt')
                    with open(ligand_out_path, 'r') as f:
                        ligand_out_data = f.read()
                    molecule = Molecule(smiles=None, input_format='pdbqt', input_data=ligand_out_data)
                    mol2_data = molecule.mol2
                    qvina_output_dict[str(id)] = mol2_data

            collection.update_one({"_id": doc["_id"]}, {"$set": {'qvina_output': qvina_output_dict}})

        # Ensure temp files are removed
        if os.path.exists(temp_path):
            os.remove(temp_path)
        if os.path.exists(ligand_path):
            os.remove(ligand_path)
        if os.path.exists(ligand_out_path):
            os.remove(ligand_out_path)
    
    def destroy(self):
        """
        Drops the database associated with this instance of the Database class.
        """
        self.client.drop_database(self.name)

    def collections_list(self):
        """
        Lists the names of all collections in the database associated with this instance of the Database class.
        """
        print('\nCollections in {}:\n'.format(self.name))
        allcollections = self.DB.list_collection_names()
        for x in allcollections:
            print(x)
        print('\n')

    def get_collection(self, COLLECTION):
        """
        Retrieves all documents from a specified MongoDB collection and converts them to a pandas DataFrame.

        Parameters:
            COLLECTION: str, name of the MongoDB collection

        Returns:
            col_df: pandas DataFrame containing all documents in the specified collection
        """
        collection = self.DB[COLLECTION]
        all = collection.find({}, {'_id': False})
        col_df = pd.DataFrame(all)
        print('\nTotal number of documents in ({collection}): {amount}'.format(
            collection=COLLECTION, amount=collection.count_documents({})) )
        return col_df

    def rm_field(self, collection_name, field_name):
        """
        Removes a specified field from all documents in a given MongoDB collection.

        Parameters:
            collection_name: str, name of the MongoDB collection
            field_name: str, name of the field to remove

        Raises:
            ValueError: If the specified collection or field does not exist
        """

        # Check if collection exists
        if collection_name not in self.DB.list_collection_names():
            raise ValueError(f"Collection '{collection_name}' does not exist in the database.")
        
        # Check if field exists in the collection
        if field_name not in self.explore(show=False)[collection_name]:
            raise ValueError(f"Field '{field_name}' does not exist in collection '{collection_name}'.")

        # Define the update query
        update_query = {"$unset": {field_name: ""}}

        # Apply the update to all documents
        collection = self.DB[collection_name]
        collection.update_many({}, update_query)

        print(f"Field '{field_name}' has been deleted from collection '{collection_name}'.")

    def rm_collection(self, collection_name):
        """
        Drops a specified collection from the database.

        Parameters:
            collection_name: str, name of the MongoDB collection to drop

        Raises:
            OperationFailure: If an error occurs when dropping the collection
        """
        try:
            collection = self.client[self.name][collection_name]
            return collection.drop()
        except OperationFailure as err:
            logging.error(f"Error dropping collection {self.name}.{collection_name}: {err}")
            raise
    
    def retain_first_n_docs(self, collection_name, N=100):
        """
        Retains only the first N documents in a specified MongoDB collection and deletes the rest.

        Parameters:
            collection_name: str, name of the MongoDB collection
            N: int, number of documents to retain (defaulted to 100)

        Raises:
            ValueError: If the specified collection does not exist
        """

        # Check if collection exists
        if collection_name not in self.DB.list_collection_names():
            raise ValueError(f"Collection '{collection_name}' does not exist in the database.")

        collection = self.DB[collection_name]

        # Get the _id of the N-th document
        nth_doc_cursor = collection.find().skip(N-1).limit(1)
        nth_doc_list = list(nth_doc_cursor)
        if len(nth_doc_list) == 0:
            print(f"The collection '{collection_name}' has fewer than {N} documents.")
            return

        nth_id = nth_doc_list[0]["_id"]
        
        # Delete all documents with _id greater than nth_id
        result = collection.delete_many({"_id": {"$gt": nth_id}})
        
        print(f"Deleted {result.deleted_count} documents. Only the first {N} documents are retained in collection '{collection_name}'.")

