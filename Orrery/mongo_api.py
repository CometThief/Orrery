import pymongo as pm
from pymongo import MongoClient
from datetime import date
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDistGeom as molDG
from tqdm.auto import tqdm
import numpy as np
import time
from collections import defaultdict
from Orrery.classes import Molecule, Database
import gzip

#DOMAIN = '172.20.0.2'
PORT = 27017

client = MongoClient(
    host = 'mongo_db',
    serverSelectionTimeoutMS = 3000, # 3 second timeout
    username = "unam",
    password = "12345",
)

#print(client.server_info())

def see(DB = None, COLLECTION = None):
    
    """
    Displays either the databases, collections within a database, or documents within a collection.

    Args:
        DB (str, optional): The name of the database. If None, list of databases will be printed.
        COLLECTION (str, optional): The name of the collection. If None, list of collections in the DB will be printed.
                                    Otherwise the contents of the table will be printed and returned as a df.

    Returns:
        pd.DataFrame: DataFrame representation of documents if DB and COLLECTION are both specified, else None.
    """
    
    col_df = None
    if DB != None:
        dbname = DB
        DB = client[DB]

    # table contents request 
    if (DB != None) and (COLLECTION != None):
        collection = DB[COLLECTION]
        all = collection.find({}, {'_id': False})
        col_df = pd.DataFrame(all)
        print('\nTotal number of documents in ({collection}): {amount}'.format(
            collection=COLLECTION, amount=collection.count_documents({})) )
        print(col_df)
        return col_df

    # list of tables within db request
    elif (DB != None) and (COLLECTION == None):
        print('\nCollections in {}:\n'.format(dbname))
        allcollections = DB.list_collection_names()
        for x in allcollections:
            print(x)
        print('\n')

    # list of dbs request
    else:
        print('\n')
        for db in client.list_databases():
            print(str(db)[1:-1])
        print('\n')

def send_zinc_tranche(tranche, DB_NAME='universe', COLLECTION = 'zinc20'):

    """
    This function takes tranche files from zinc20 and parses them then sends them to Mongo.
    Any new filetypes must be implemented.

    Default schema used:

    document = {
                            'smiles_chain' : smiles_chain,
                            'zinc_id' : zinc_id,
                            'tranche' : tranche,
                            'atomspresent' : atomspresent
                        }

    Filetypes accepted so far: 
    .smi
    """
    
    batch = []
    working_db = Database(DB_NAME)
    
    # parsing many-molecule data from .smi files
    if tranche.endswith('smi'):
        # organize the molecules in a single tranche
        with open(tranche, 'rt') as file:
            tranche_string = file.read()
        lines = tranche_string.strip().split('\n')

        #
        for line in lines[1:]:
            smiles, zinc_id = line.split()
            mol = Molecule(smiles, zinc_id)
            document = {
                            'smiles' : smiles,
                            'mol2' : mol.mol2,
                            'zinc_id' : zinc_id
                        }
            batch.append(document)
    else:
        __, file_extension = os.path.splitext(tranche)
        print(f'Unsupported file type: {file_extension}')
        pass

    working_db.insert(COLLECTION, batch)

def subsearch(DB_NAME, COLLECTION, substructs = pd.DataFrame(
        {
            'name': ['carbonyl'],
            'smiles': ['C=O']
        }
    ) , newname = None, showall=False):

    """
    Searches for molecules matching a substructure pattern in a MongoDB collection.

    Args:
        DB_NAME (str): The name of the database.
        COLLECTION (str): The name of the collection.
        substructs (pd.DataFrame, optional): DataFrame with name and SMILES string of substructures. Defaults to carbonyl.
        newname (str, optional): Name of the new collection to insert matched molecules. Defaults to 'subsearch' with current timestamp.
        showall (bool, optional): If True, prints all matched molecules. Defaults to False.
    """

    if not newname:
        newname = 'subsearch ' + time.strftime("%Y%m%d-%H%M%S")

    working_db = Database(DB_NAME)
    substructs['rdkit molecule'] = substructs['smiles'].apply(Chem.MolFromSmiles)
    col = see(DB_NAME, COLLECTION)
    matches = []

    for index, row in tqdm(col.iterrows(), total=col.shape[0]):
        mol = Chem.MolFromSmiles(row['smiles'])
        match = False
        substructmatches = list()
        for _, substruct in substructs.iterrows():
            if mol:
                if mol.HasSubstructMatch(substruct['rdkit molecule']):
                    substructmatches.append(substruct['name'])
                    match = True
        if match == True:
            matches.append(
                    {
                        'smiles': row['smiles'],
                        'substructure match': substructmatches
                    }
                )
        #if not match:
            #no_matches.append(index)

    if matches and showall==True:
        with pd.option_context('display.max_rows', None,
                       'display.max_columns', None,
                       'display.precision', 3,
                       ):
            print(matches)
    elif matches:
        working_db.insert(newname, matches)
        matches=pd.DataFrame(matches)
        #no_matches=pd.DataFrame(no_matches)
        print('{0} molecules match the filters'.format(len(matches)))
        #print('{0} molecules dont'.format(len(no_matches)))
    else:
        print('\nNo values matched the filters\n')

# UNDER MAINTENANCE
def sizefilter(DB_NAME, COLLECTION, min = 3, max = 7, newname = None):

    """
    Filters molecules based on a size criteria from a MongoDB collection.

    Args:
        DB_NAME (str): The name of the database.
        COLLECTION (str): The name of the collection.
        min (int, optional): Minimum size criteria. Defaults to 3.
        max (int, optional): Maximum size criteria. Defaults to 7.
        newname (str, optional): Name of the new collection to insert filtered molecules. Defaults to formatted string with COLLECTION, min, max, and date.
    """

    db = client[DB_NAME]
    collection = db[COLLECTION]
    filtered = []
    if not newname:
        newname = 'Size filter from ({collection}): > {min} ; < {max}  {date}'.format(
            collection=COLLECTION, min=min, max=max, date=time.strftime("%Y%m%d-%H%M%S"))

    filtercursor = collection.find( {'$and': [
        { 'maxdist': { '$gt': min } },
        { 'maxdist': { '$lt': max } } 
        ] } )
    
    for x in filtercursor:
        filtered.append(x)
    
    insert(DB_NAME, newname, filtered)
    print('Total number of molecules that fit the size criteria: ', len(filtered))


def atomfilter(DB_NAME, COLLECTION, contains = [7, 17], newname = None):
    #, exclude=range(min,max)

    """
    Filters molecules based on the presence of certain atoms from a MongoDB collection.

    Args:
        DB_NAME (str): The name of the database.
        COLLECTION (str): The name of the collection.
        contains (List[int]): List of atomic numbers that must be present in a molecule.
        newname (str, optional): Name of the new collection to insert filtered molecules. Defaults to formatted string with contains, COLLECTION, and date.
    """
    working_db = Database(DB_NAME)
    collection = db[COLLECTION]
    if not newname:
        newname = '{contains} Atom filter from ({COLLECTION})  {DATE}'.format(
            contains=contains, COLLECTION=COLLECTION, DATE=time.strftime("%Y%m%d-%H%M%S")
        )

    print('\nFiltering from {DB_NAME} - {COLLECTION}:'.format(
        DB_NAME=DB_NAME, COLLECTION=COLLECTION
    ))

    allvalues = collection.find()
    results = []

    for value in allvalues:
        atoms = value['atomspresent']
        check = list(set(atoms) - set(contains))
        
        if not check:
            results.append(value)

    if results:
        print('{length} values met the criteria'.format(length = len(results)))
        insert(DB_NAME, newname, results)
    else: 
        print('\nNo values met the criteria')

def to_gaussinput(DB_NAME, COLLECTION):

    """
    Converts a MongoDB collection of molecules into Gaussian input file format.

    Args:
        DB_NAME (str): The name of the database.
        COLLECTION (str): The name of the collection.
    """

    db = client[DB_NAME]
    collection = db[COLLECTION]
    allvalues = pd.DataFrame(collection.find())
    directory = './gauss_files/input/input.smi'
    content = '\n'.join((allvalues['smiles'].tolist()))

    with open(directory, 'w') as f:
        f.write(content)

DB = 'dna'
collname = 'smiles ' + str(date.today())

