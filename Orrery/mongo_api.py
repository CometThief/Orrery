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
from Orrery.classes import mol2
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

def insert(DB, COLLECTION, documents):

    """
    Inserts multiple documents into a specified database / collection in MongoDB.

    Args:
        DB (str): The name of the database.
        COLLECTION (str): The name of the collection.
        documents (List[dict]): A list of documents to insert.

    Returns:
        pymongo.results.InsertManyResult: Object with details of the insert operation.
    """

    DB = client[DB]
    COLLECTION = DB[COLLECTION]
    inserted = COLLECTION.insert_many(documents)
    return inserted

def destroy(DB):

    """
    Drops a specified database from MongoDB.

    Args:
        DB (str): The name of the database to drop.
    """

    client.drop_database(DB)

def rm_collection(DB, COLLECTION):

    """
    Drops a specified collection from a specified database in MongoDB.

    Args:
        DB (str): The name of the database.
        COLLECTION (str): The name of the collection to drop.

    Returns:
        bool: True if collection dropped successfully, else False.
    """

    DB = client[DB]
    COLLECTION = DB[COLLECTION]
    deleted = COLLECTION.drop()
    return deleted

def see(DB = None, COLLECTION = None):
    
    """
    Displays either the databases, collections within a database, or documents within a collection.

    Args:
        DB (str, optional): The name of the database. If None, list of databases will be printed.
        COLLECTION (str, optional): The name of the collection. If None, list of collections in the DB will be printed.

    Returns:
        pd.DataFrame: DataFrame representation of documents if DB and COLLECTION are both specified, else None.
    """
    
    col_df = None

    if DB != None:
        dbname = DB
        DB = client[DB]

    if (DB != None) and (COLLECTION != None):
        collection = DB[COLLECTION]
        all = collection.find({}, {'_id': False})
        col_df = pd.DataFrame(all)

        print('\nTotal number of documents in ({collection}): {amount}'.format(
            collection=COLLECTION, amount=collection.count_documents({})) )

    elif (DB != None) and (COLLECTION == None):
        print('\nCollections in {}:\n'.format(dbname))
        allcollections = DB.list_collection_names()
        for x in allcollections:
            print(x)
        print('\n')

    else:
        print('\n')
        for db in client.list_databases():
            print(str(db)[1:-1])
        print('\n')

    if col_df: print(col_df)
    return col_df

def parse_tranche(tranche, DB_NAME='universe', COLLECTION = 'smiles ' + time.strftime("%Y%m%d-%H%M%S")):

    """
    This function takes tranche files from zinc20 and parses them. Any new filetypes must be implemented.

    Filetypes accepted so far: 
    .mol2.gz
    """
    batch = {}
    
    # parsing many-molecule data from .mol2.gz files
    if tranche.endswith('mol2.gz'):
        with gzip.open(tranche, 'rt') as file:
            mol_string = file.read()
        molecules = []
        mol_data = mol_string.split('@<TRIPOS>MOLECULE')
        for data in mol_data[1:]:
            molecule = mol2('@<TRIPOS>MOLECULE' + data)
            molecules.append(molecule)
        for mol in molecules:
            print(mol)
            print(mol.sections)
            print(mol.mol_string)
            input('heee')
            #batch.update()
    else:
        __, file_extension = os.path.splitext(tranche)
        print('Unsupported file type: {file_extension}')
        pass

    # here one might add other parsing conditions for other file extensions
    print(molecules)


    
    '''
    else:
        print(f"The smiles directory is empty. Attemps left: '{maxemptyiterations - iterations}'")
        iterations += 1
        time.sleep(5)
    else:
        for i in allfiles:
            tranche = i[-10:]
            with open(i, 'r') as file:
                for e in file:
                    #this should be generalized with regex to be able to handle more file types (?)
                    #only works with zinc15 smiles files as is
                    if not e.startswith('smiles'):
                        smiles_chain = e.split()[0]
                        zinc_id = e.split()[1]
                        mol = Chem.MolFromSmiles(smiles_chain)

                        #distance matrix stuff
                        distmatrix = molDG.GetMoleculeBoundsMatrix(mol)
                        maxdist = np.max(distmatrix)

                        #atoms present in molecule by atomic number
                        atomic_count = defaultdict(lambda : 0)
                        for atom in mol.GetAtoms():
                            atomic_count[atom.GetAtomicNum()] += 1
                        atomspresent = [i for i in atomic_count if atomic_count[i]!=atomic_count.default_factory()]

                        document = {
                            'smiles_chain' : smiles_chain,
                            'zinc_id' : zinc_id,
                            'tranche' : tranche,
                            'maxdist' : maxdist,
                            'atomspresent' : atomspresent
                        }
                        batch.append(document)
                        counter += 1
            iterations = 0
            #print('Finished tranche: {ctranche}'.format(ctranche=i))
            os.remove(i)
            to_remove.append(i)

        if to_remove:
            for i in to_remove:
                allfiles.remove(i)
        to_remove = list()

        # decided to place this inside the loop so many mini batches are inserted instead of one large one, this could be a mistake       
        if batch:
            insert(DB_NAME, COLLECTION, batch)
            batch = list()
        '''
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

    substructs['rdkit molecule'] = substructs['smiles'].apply(Chem.MolFromSmiles)

    #mols = list(substructs['rdkit molecule'])
    #names = list(substructs['name'])
    #output image with substructures?
    #Chem.Draw.MolsToGridImage(
    # mols,
    # legends=name
    # molsPerRow=5
    # )

    col = see(DB_NAME, COLLECTION)
    matches = []
    #no_matches = []

    for index, row in tqdm(col.iterrows(), total=col.shape[0]):
        mol = Chem.MolFromSmiles(row['smiles_chain'])
        match = False
        substructmatches = list()
        for _, substruct in substructs.iterrows():
            if mol.HasSubstructMatch(substruct['rdkit molecule']):
                substructmatches.append(substruct['name'])
                '''
                matches.append(
                    {
                        'zinc_id': row['zinc_id'],
                        'smiles_chain': row['smiles_chain'],
                        'substructure match': substruct['name']
                    }
                )
                '''
                match = True
        if match == True:
            matches.append(
                    {
                        'zinc_id': row['zinc_id'],
                        'smiles_chain': row['smiles_chain'],
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
        insert(DB_NAME, newname, matches)
        matches=pd.DataFrame(matches)
        #no_matches=pd.DataFrame(no_matches)
        print('{0} molecules match the filters'.format(len(matches)))
        #print('{0} molecules dont'.format(len(no_matches)))
    else:
        print('\nNo values matched the filters\n')


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


def atomfilter(DB_NAME, COLLECTION, contains, newname = None):
    #, exclude=range(min,max)

    """
    Filters molecules based on the presence of certain atoms from a MongoDB collection.

    Args:
        DB_NAME (str): The name of the database.
        COLLECTION (str): The name of the collection.
        contains (List[int]): List of atomic numbers that must be present in a molecule.
        newname (str, optional): Name of the new collection to insert filtered molecules. Defaults to formatted string with contains, COLLECTION, and date.
    """

    db = client[DB_NAME]
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

