#This file contains the definition for the Molecule class which is used to represent and process molecular data

import pandas as pd
from pymongo import MongoClient

class mol2:
    def __init__(self, mol_string):
        self.mol_string = mol_string
        self.sections = self._parse_mol2()

    def _parse_mol2(self):
        sections = self.mol_string.split('@<TRIPOS>')
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