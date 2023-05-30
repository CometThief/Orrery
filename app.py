from Orrery.chemcollector import fetch_zinc, fetch_chembl, fetch_REAL
from Orrery.mongo_api import see, destroy
#import smiles_gaussian
import pandas as pd
import time

def fetch_all():

    fetch_zinc()
    #fetch_chembl()
    #fetch_REAL()

def Abaddon(database = 'universe'):
    """
    The Hebrew term Abaddon (Hebrew: אֲבַדּוֹן ’Ăḇaddōn, meaning "destruction", "doom"), 
    and its Greek equivalent Apollyon (Koinē Greek: Ἀπολλύων, Apollúōn meaning "Destroyer") 
    appear in the Bible as both a place of destruction and an angel of the abyss. In the 
    Hebrew Bible, abaddon is used with reference to a bottomless pit, often appearing alongside 
    the place Sheol (שְׁאוֹל Šəʾōl), meaning the resting place of dead peoples.

            -Wikipedia
    """
    destroy(database)

def main():
    see()
    fetch_all()
    input('here?')

if __name__ == "__main__":
    #fetch_all()
    main()

# download from zinc20

#mapi.rm_db('universe')
#mapi.rm_collection('universe', 'chembl')
#col = mapi.visualize_contents('universe', 'chembl')
#cc.fetch_chembl(limit=5)
#print(col.columns)
#cc.fetch_all()
#mapi.visualize_contents()
#print(col0)
#see('universe', 'REAL')
#col1 = mapi.visualize_contents('universe', 'zinc20')
#print(col1)
#mapi.rm_db('universe')
#col = mapi.visualize_contents('universe', 'chembl')
#print(col)
#cc.fetch_all()
#cc.fetch_zinc()
#cc.fetch_chembl_1()
#print('continuing')

# grab downloaded files from zinc15 and send them to db
#mapi.rm_db('Zinc15')
#time.sleep(10)
#mapi.grab_smiles('universe', 'Zinc20')

# write smiles file from specific collection
#mapi.to_gaussinput('Zinc15', 'subsearch 20221213-183536')

# run transformation from 
#smiles_gaussian.gauss_gen()

# delete a DB or collection
#mapi.rm_db('Zinc15')
#mapi.rm_collection('Zinc15', 'Size filter from (Zinc15 Universe): > 3 ; < 7  20230110-014558')
#mapi.rm_collection('Zinc15', 'Test03')

# query - not working
#mapi.query('dna','smiles 2022-11-04')

# size filter
#mapi.sizefilter('Zinc15', 'Zinc15 Universe', min=3, max=7, newname = 'Test01')

# atom filtering
#mapi.atomfilter('Zinc15', 'Test01', contains=(6,7,8), newname = 'Test02')

# substructure search
'''
mapi.subsearch('Zinc15', 'Test02', substructs = pd.DataFrame(
        {
            'name': ['carbonyl', 'propyl'],
            'smiles': ['C=O','CCC']
        }),
        newname = 'Test03')
'''
# visualization
#collection = mapi.visualize_contents()
#print(collection)
#print(type(collection))


#{'smiles': 'C[C@@H](CC(=O)O)NC(=O)OCC1C2=CC=CC=C2C2=CC=CC=C21', 'idnumber': 'Z1982491139', 'MW': '325.364', 'HAC': '24', 'sLogP': '3.388', 'HBA': '3', 'HBD': '2', 'RotBonds': '5', 'FSP3': '0.263', 'TPSA': '75.630', 'QED': '0.883', 'PAINS': '', 'BRENK': '', 'NIH': '', 'ZINC': '', 'LILLY': '', 'lead-like': 'True', '350/3_lead-like': '', 'fragments': '', 'strict_fragments': '', 'PPI_modulators': '', 'natural_product-like': '', 'Type': 'S', 'InChiKey': 'LYMLSPRRJWJJQD-LBPRGKRZSA-N'}

