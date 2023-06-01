from Orrery.chemcollector import fetch_zinc, fetch_chembl, fetch_REAL
from Orrery.mongo_api import see, subsearch
from Orrery.classes import Database
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
    ''''''
    ##### some tests

    ### cloning a db locally
    #fetch_zinc()
    #fetch_chembl()
    #fetch_REAL()
    #fetch_all()

    ### see() usage
    #see() # list all databases
    #see('universe') # list all collections in 'universe' db
    #see('universe', 'REAL') # visualize a df containing the collection 'REAL' which is inside the 'universe' db
    #realdf = see('universe', 'REAL') # save the df to be used later

    ### subsearching - searches for C=O by default

    # dict struct for subsearching:
    '''
    {
            'name': ['carbonyl'],
            'smiles': ['C=O']
        }
    '''

    #subsearch('universe', 'chembl')
    #subsearch('universe', 'REAL')
    #subsearch('universe', 'zinc20')

    ### destroy a db
    #Abaddon()
    #input('see')

    
if __name__ == "__main__":
    main()
