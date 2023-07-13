from Orrery.chemcollector import fetch_zinc, fetch_chembl, fetch_REAL
from Orrery.mongo_api import see, subsearch
from Orrery.classes import Database, Molecule
#import smiles_gaussian
import pandas as pd
import time
from rdkit import Chem

def fetch_all():

    fetch_zinc()
    fetch_chembl()
    fetch_REAL()

def Abaddon(database = 'universe'):
    """
    The Hebrew term Abaddon (Hebrew: אֲבַדּוֹן ’Ăḇaddōn, meaning "destruction", "doom"), 
    and its Greek equivalent Apollyon (Koinē Greek: Ἀπολλύων, Apollúōn meaning "Destroyer") 
    appear in the Bible as both a place of destruction and an angel of the abyss. In the 
    Hebrew Bible, abaddon is used with reference to a bottomless pit, often appearing alongside 
    the place Sheol (שְׁאוֹל Šəʾōl), meaning the resting place of dead peoples.

            -Wikipedia
    """
    goodbye = Database(database)
    goodbye.destroy()

def main():

    ################################
    ##### cloning a db locally #####
    ################################

    
    '''
    # clean existing db
    Abaddon('universe')
    # individual fetching - would refetching be equivalent to updating?
    fetch_chembl()
    fetch_chembl()
    fetch_REAL()
    
    # grouped fetching
    fetch_all()
    '''


    #######################
    ##### see() usage #####
    #######################


    #see() # list all databases
    #see('universe') # list all collections in 'universe' db
    #see('universe', 'chembl') # visualize a df containing the collection 'REAL' which is inside the 'universe' db
    #realdf = see('universe', 'chembl') # save the df to be used later
    #print(realdf.head(20))


    
    ###########################
    ##### common querying #####
    ###########################


    '''
    # this is an example of how to implement basic querying, it needs to be heavily refined
    # this syntax is specific to mongo
    wdb = Database('universe')
    # .explore() prints out a bunch of information about the table being accessed
    wdb.explore()
    # convert the column from string to double
    wdb.convert_collection_column('chembl', 'alogp', 'double')
    # filter by all values of alogp greater than 2.1
    results = wdb.query('chembl', {'alogp':{'$gt':2.1}}, store_results = False)
    for x in results:
        print(x)
    '''
    
    
    ####################################
    ##### minimal qvina processing #####
    ####################################


    # clean the db first  then fetch 1 single item
    #Abaddon('universe')
    #fetch_chembl(limit=1, max_pages=1)
    
    # generate a mol2 column that is ph-adjusted by default
    working_db = Database('universe')
    print(f'---Converting the smiles column from the chembl collection to a dict of pH-adjusted mol2 files---')
    working_db.convert_smiles_to_mol2('chembl')
    print(f'---Converting the pH_n columns to qvina outputs---')
    working_db.basic_qvina_analysis_timed('chembl')

    # quick and dirty dataframe testing - not recommended for large datasets
    chembl_df = see('universe','chembl')
    single_item = chembl_df['qvina_output'][0]
    for pH, qvina_file in single_item.items():
        print(f'pH {pH}:\n{qvina_file}\n\n')
    
if __name__ == "__main__":
    main()
    
