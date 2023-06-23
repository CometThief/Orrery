#This file contains the definition for the Molecule class which is used to represent and process molecular data

import pandas as pd
from pymongo import MongoClient, errors
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem
from openbabel import openbabel
import pybel
import logging
from datetime import datetime
import math

def compute_scaffold(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smiles = Chem.MolToSmiles(scaffold)
        return scaffold_smiles
    except Exception as e:
        print(f"Failed to compute scaffold for SMILES {smiles}: {e}")
        return math.nan

class Molecule:
    def __init__(self, smiles, id=None):
        self.smiles = smiles
        self.mol_id = id
        self.mol2 = self._get_mol2()
        self.sections = self._parse_mol2()
        self.dict = self.to_dict()

    def _get_mol2(self):
        mol = Chem.MolFromSmiles(self.smiles)
        mol = Chem.AddHs(mol) # Add Hs
        AllChem.EmbedMolecule(mol) # Generate 3D coordinates
        AllChem.UFFOptimizeMolecule(mol) # universal force field optimization
        obmol = openbabel.OBMol()
        converter = openbabel.OBConversion()
        converter.SetInAndOutFormats("mol", "mol2")
        converter.ReadString(obmol, Chem.MolToMolBlock(mol))
        if self.mol_id : obmol.SetTitle(self.mol_id) # set optional id (replaced with asterisks if 'None')
        mol2_string = converter.WriteString(obmol) # Get mol2 string
        return mol2_string

    def _parse_mol2(self):
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

    def to_dict(self):
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
        try:
            self.client = MongoClient(
                   host = host,
                port = port,
                serverSelectionTimeoutMS = 3000, # 3 second timeout
                username = username,
                password = password,
            )
            self.client.server_info()  # Trigger ServerSelectionTimeoutError if cannot connect to MongoDB server
        except ServerSelectionTimeoutError as err:
            logging.error(f"Could not connect to MongoDB: {err}")
            raise
        self.DB = self.client[self.name]

    def insert(self, collection_name, documents):
        collection = self.client[self.name][collection_name]
        collection.create_index("smiles", unique=True)
        try:
            collection.insert_many(documents, ordered=False)
        except errors.BulkWriteError as e:
            #print(e.details['writeErrors'])
            ''

    def explore(self, show=True):
        collections_info = {}
        collections = self.DB.list_collection_names()
        print("\n--- Collections Information ---")
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

    def query(self, collection_name, query, store_results=True):
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
        '''
        Method to convert an entire column in a mongo collection a certain data type.

        To convert to floats, pass the argument 'double' instead of 'float'.

        If conversion fails, original entry will be kept by default.

        To convert a float stored as a str to int, it's necessary to convert to 'double' first then to 'int'.
        '''
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
        collection = self.DB[collection]
        documents = collection.find()
        
        for document in documents:
            smiles = document["smiles"]
            scaffold_smiles = compute_scaffold(smiles)
            collection.update_one({"_id": document["_id"]}, {"$set": {"scaffold": scaffold_smiles}})
        
        if index: collection.create_index("scaffold")

    def scaffold_search(self, collection_name, smiles, store_results=True):
        
        # Make sure the 'scaffold' column exists
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

    def destroy(self):
        self.client.drop_database(self.name)

    def collections_list(self):      
        print('\nCollections in {}:\n'.format(self.name))
        allcollections = self.DB.list_collection_names()
        for x in allcollections:
            print(x)
        print('\n')

    # this method returns the contents of the collection as a df
    def get_collection(self, COLLECTION):
        collection = self.DB[COLLECTION]
        all = collection.find({}, {'_id': False})
        col_df = pd.DataFrame(all)
        print('\nTotal number of documents in ({collection}): {amount}'.format(
            collection=COLLECTION, amount=collection.count_documents({})) )
        return col_df

    def rm_collection(self, collection_name):
        try:
            collection = self.client[self.name][collection_name]
            return collection.drop()
        except OperationFailure as err:
            logging.error(f"Error dropping collection {self.name}.{collection_name}: {err}")
            raise



'''
class Database:
    def __init__(self, name, client):
        self.name = name
        self.client = client

    def get_collection(self, collection_name):
        return Collection(self.name, collection_name, self.client)

class Collection:
    def __init__(self, db_name, collection_name, client):
        self.db_name = db_name
        self.collection_name = collection_name
        self.client = client

    def insert(self, document):
        print(f"Inserting {document} into {self.collection_name} in {self.db_name} database.")

    def fetch(self, query={}):
        print(f"Fetching documents from {self.collection_name} in {self.db_name} database with query {query}.")

    def remove(self, query):
        print(f"Removing documents from {self.collection_name} in {self.db_name} database with query {query}.")

    def check_duplicate(self, key):
        print(f"Checking for duplicates in {self.collection_name} in {self.db_name} database with key {key}.")
'''

#working_db = Database('hola')
#working_db.substruct()