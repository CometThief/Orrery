"""
chemcollector.py

This module contains functions for fetching and parsing chemical data from the Zinc20, ChEMBL, and REAL chemical databases.
It facilitates web scraping, data extraction and conversion, as well as data loading into a MongoDB collection for subsequent use.
The chemcollector module relies on the requests, pandas, tqdm, xml, pickle, and rdkit libraries to perform these operations.
"""

from re import L
import requests
import math
import os
from tqdm import tqdm
from . import mongo_api as mapi
import multiprocessing
import requests
import xml.etree.ElementTree as ET
import pandas as pd
import re
import bz2
import shutil
from rdkit import Chem
import pickle
from .classes import Database
import tarfile


def fetch_zinc(rep="3D", since="", db_r="", format="smi",
               using="uri", DB_NAME = 'universe', COLLECTION='zinc20',
               all_3d_url = "https://zinc20.docking.org/tranches/all3D.json",
               download_url = "https://zinc20.docking.org/tranches/download"):

    """
    Fetches the entire Zinc20 database in a specified format and loads the data into a MongoDB collection.

    Parameters
    ----------
    rep : str, optional
        Representation of the molecules, by default "3D".
    since : str, optional
        To fetch data since a specific date (in YYYY-MM-DD format), by default "".
    db_r : str, optional
        To limit data by database release, by default "".
    format : str, optional
        The format of the molecular data to be fetched, by default "mol2.gz".
    using : str, optional
        Method of fetching data, by default "uri".
    DB_NAME : str, optional
        Name of the MongoDB database to load data into, by default 'universe'.
    COLLECTION : str, optional
        Name of the MongoDB collection to load data into, by default 'zinc20'.

    Example
    -------
    fetch_zinc(DB_NAME="mydatabase", COLLECTION="mycollection")

    Notes
    -----
    This function fetches data downloading individual Zinc tranches. To keep track of the progress, a list of finished tranches is maintained.
    """

    r = requests.get(all_3d_url)
    tranches = r.json()

    # this is the format used by zinc since at least zinc15
    data = {
        'representation': rep,
        'tranches': ' '.join(x['name'] for x in tranches),
        'format': format,
        'using': using
    }

    print('Requesting zinc20...')
    r = requests.post(download_url, data=data, stream=True)

    # Create the dir for the tranches
    WORKING_DIR = './tranches/'
    os.makedirs(WORKING_DIR, exist_ok=True)

    #the '0_' is so the 2 files line up at the top of the 

    with open(WORKING_DIR + '0_alltranches', 'w') as a:
        total = r.text
        a.write(total)
        total = len(str(total).split('\n'))
    
    try:
        with open(WORKING_DIR + '0_finished_tranches', 'r') as all_file:
            done_list = set(all_file.read().split())
    except FileNotFoundError:
        done_list = set()

    print('Downloading Zinc20 tranches...')
    for x in tqdm(r.iter_lines(), total=total):

        if (str(x)) not in done_list:
            #getting each tranche and writing the gz file
            resp = requests.get(x)
            # for some reason some links are dead. probably internal zinc inconsistencies. this if(resp.ok) sidesteps the problem.
            if(resp.ok):
                path = x.decode('utf-8')
                current_file_name = os.path.join(WORKING_DIR, os.path.basename(path))
                with open (current_file_name, 'wb+') as f:
                    f.write(resp.content)

                #writes to a file with a list of correctly downloaded tranches
                with open(WORKING_DIR + '0_finished_tranches', 'a') as f:
                    f.write(str(x) + '\n')
                
                mapi.send_zinc_tranche(current_file_name, DB_NAME, COLLECTION)

def fetch_chembl(DB_NAME = 'universe', COLLECTION = 'chembl', limit = 1000, rdkit_obj = True, max_pages = None):
    """
    Fetches a specified number of molecules from the ChEMBL database and loads the data into a MongoDB collection.
    Fetches a certain amount of molecules per a page for a certain amount of pages.

    Parameters
    ----------
    DB_NAME : str, optional
        Name of the MongoDB database to load data into, by default 'universe'.
    COLLECTION : str, optional
        Name of the MongoDB collection to load data into, by default 'chembl'.
    limit : int, optional
        Number of molecules to fetch per page, by default 1000.
    rdkit_obj : bool, optional
        Whether to convert molecular structures into RDKit object and include them in the data, by default True.
    max_pages : int, optional
        Maximum number of pages to scrape. If not provided or None, all pages are scraped.

    Example
    -------
    fetch_chembl(DB_NAME="mydatabase", COLLECTION="mycollection", limit=500, max_pages=5)

    """

    # TO DO: stop/restart where left off; dupe protection

    # Define the base URL and the initial URL
    base_url = 'https://www.ebi.ac.uk'
    limit = str(limit)
    url = base_url + '/chembl/api/data/molecule?limit=' + limit
    chembl_db = Database(DB_NAME)
    batch = list()
    print('Processing chembl database...')
    counter = 1
    while url:
        if max_pages is not None and counter > max_pages:
            break
        response = requests.get(url)
        root = ET.fromstring(response.content)
        print(f'\rScraping page "{counter}" of REAL...', end='', flush=True)
        xml_str = response.content.decode('utf-8')
        batch = parse_xml_str(xml_str, rdkit_obj)
        counter+=1
        url = base_url + root.find('page_meta').find('next').text if root.find('page_meta').find('next') is not None else None
        chembl_db.insert(COLLECTION, batch)
        batch = list()

def fetch_REAL(username = 'samuelmar@gmail.com', password = 'TCa!@wjrrH9F4Jv', directory = './temp/'):
    
    """"
    Downloads the Enamine REAL database, unzips the files, and loads the data into a MongoDB collection.

    Parameters
    ----------
    username : str, optional
        Username for the Enamine account, by default 'samuelmar@gmail.com'.
    password : str, optional
        Password for the Enamine account, by default 'TCa!@wjrrH9F4Jv'.
    directory : str, optional
        Directory to store the downloaded files, by default './temp/'.

    Example
    -------
    fetch_REAL(username='yourusername', password='yourpassword', directory='./yourdirectory/')

    Notes
    -----
    This function downloads large bz2 files. Progress is displayed in the terminal.
    """

    if not os.path.exists(directory):
        os.makedirs(directory)

    login_url = 'https://enamine.net/component/users/?task=user.login'
    download_urls = (
        'https://ftp.enamine.net/download/REAL/Enamine_REAL_HAC_24_394M_CXSMILES.cxsmiles.bz2',
        'https://ftp.enamine.net/download/REAL/Enamine_REAL_HAC_22_23_471M_CXSMILES.cxsmiles.bz2',
        'https://ftp.enamine.net/download/REAL/Enamine_REAL_HAC_24_394M_CXSMILES.cxsmiles.bz2',
        'https://ftp.enamine.net/download/REAL/Enamine_REAL_HAC_25_557M_CXSMILES.cxsmiles.bz2',
        'https://ftp.enamine.net/download/REAL/Enamine_REAL_HAC_26_833M_Part_1_CXSMILES.cxsmiles.bz2',
        'https://ftp.enamine.net/download/REAL/Enamine_REAL_HAC_26_833M_Part_2_CXSMILES.cxsmiles.bz2',
        'https://ftp.enamine.net/download/REAL/Enamine_REAL_HAC_27_1.1B_Part_1_CXSMILES.cxsmiles.bz2',
        'https://ftp.enamine.net/download/REAL/Enamine_REAL_HAC_27_1.1B_Part_2_CXSMILES.cxsmiles.bz2',
        'https://ftp.enamine.net/download/REAL/Enamine_REAL_HAC_28_1.2B_Part_1_CXSMILES.cxsmiles.bz2',
        'https://ftp.enamine.net/download/REAL/Enamine_REAL_HAC_28_1.2B_Part_2_CXSMILES.cxsmiles.bz2',
        'https://ftp.enamine.net/download/REAL/Enamine_REAL_HAC_29_38_988M_Part_1_CXSMILES.cxsmiles.bz2',
        'https://ftp.enamine.net/download/REAL/Enamine_REAL_HAC_29_38_988M_Part_2_CXSMILES.cxsmiles.bz2'
    )
    
    # set up the session
    with requests.Session() as session:
        
        # perform login
        login_payload = {
            'username': username,
            'password': password,
            'option': 'com_users',
            'task': 'user.login',
            'return': 'aHR0cHM6Ly9lbmFtaW5lLm5ldC9jb21wb3VuZC1jb2xsZWN0aW9ucy9yZWFsLWNvbXBvdW5kcy9yZWFsLWRhdGFiYXNlbw==',
            'b9d45469bc1110a7adec048f14c2bd90': '1'
        }
        # just a cute identifier of what email is being used to log in to REAL
        username, domain = username.split('@')
        masked_username = '*' * len(username[:-3]) + username[-3:]
        print(f'Logging into REAL using: {masked_username}@{domain}')
        response = session.post(login_url, data=login_payload)
        
        # check if login was successful
        if response.status_code != 200:
            print('Login to REAL failed')
            exit()
        print('Login to REAL successful')

        for url in download_urls:

            filename = re.search(r'(?<=\/)[^\/]+\.bz2$', url).group(0)
            print('Downloading file {filename}'.format(filename=filename))
            response = session.get(url, stream=True)

            # prepare the progress bar, since downloads are pretty large
            total_size_in_bytes = int(response.headers.get('content-length', 0))    
            block_size = 1024 # 1 Kibibyte
            

            # save the file to disk
            filename = directory + filename
            if os.path.exists(filename) and os.path.getsize(filename) == total_size_in_bytes:
                print(f"The file '{filename}' already exists. Skipping download.")
            else:
                progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
                with open(filename, 'wb') as file:
                    for data in response.iter_content(block_size):
                        progress_bar.update(len(data))
                        file.write(data)
                print(f'File {filename} downloaded successfully')
                progress_bar.close()
            parse_bz2(filename, directory)
            #process_REAL = multiprocessing.Process(target=parse_bz2, args=(filename,directory))
            #process_REAL.start()
               
    print('Finished cloning REAL database locally')
    os.remove(directory)

def parse_xml_str(xml_str, rdkit_obj):

    """
    This function parses an XML string into a list of molecule dictionaries. Each dictionary represents a molecule with its properties and structures.

    Args:
    xml_str (str): A string in XML format to be parsed.
    rdkit_obj (bool): A boolean indicating whether to create an RDKit object from the canonical SMILES of the molecule.

    Returns:
    list: A list of dictionaries, where each dictionary represents a molecule with its properties and structures.

    Example
    -------
    molecules = parse_xml_str(your_xml_string, rdkit_obj=True)

    """

    root = ET.fromstring(xml_str)
    # Create a list to hold the dictionaries
    molecules_list = []
    # Loop through each molecule
    for molecule in root.findall('./molecules/molecule'):
        molecule_dict = {}
        properties = molecule.find('molecule_properties')
        for prop in properties:
            if prop.text:
                molecule_dict[prop.tag] = prop.text

        structures = molecule.find('molecule_structures')
        for struct in structures:
            if struct.text:
                if struct.tag == 'canonical_smiles': # rename the smiles column to maintain consistency with other dbs
                    molecule_dict['smiles'] = struct.text
                else:
                    molecule_dict[struct.tag] = struct.text
        if rdkit_obj and 'smiles' in molecule_dict: # this condition ignores empty molecules in the chembl db
            rdkit_mol = Chem.MolFromSmiles(molecule_dict['smiles'])
            serialized_mol = pickle.dumps(rdkit_mol)
            # de-pickling:
            #rdkit_mol = pickle.loads(serialized_mol)
            molecule_dict['pickled_rdkit_object'] = serialized_mol
            molecules_list.append(molecule_dict)

    return molecules_list

def parse_bz2(filename, directory, DB_NAME = 'universe', COLLECTION = 'REAL', remove = True, max_batch_size = 1000):
    """
    This function reads a bz2 file, parses its contents, and stores the data in a MongoDB database.

    Args:
    filename (str): The path of the bz2 file.
    directory (str): The directory to store the extracted data.
    DB_NAME (str): The name of the MongoDB database, default is 'universe'.
    COLLECTION (str): The name of the MongoDB collection, default is 'REAL'.
        remove (bool): Whether to remove the bz2 file after extraction, default is True.
    max_batch_size (int): The maximum number of molecules to include in each insert operation, default is 1000.

    Returns:
    None

    Example
    -------
    parse_bz2(filename="yourfile.bz2", directory="./yourdirectory/", DB_NAME="mydatabase", COLLECTION="mycollection", remove=False, max_batch_size=500)

    Notes
    -----
    This function should probably only handle the bz2 file and the actual writing to mongo should be done elsewhere.
    This function reads a bz2 file, which typically contains a large amount of data. It also uses a batch insertion method to efficiently load data into the MongoDB database.
    """

    working_db = Database(DB_NAME)
    batch = []
    print(f"Moving '{filename}' to local database")
    with bz2.open(filename, 'rt') as file:
        keys = file.readline().strip().split("\t")
        file_size = os.path.getsize(filename)
        chunk_size = 1024 * 1024 # 1 MB
        num_chunks = file_size // chunk_size + 1
        for i in tqdm(range(num_chunks), desc=f"Processing {filename}", unit="MB"):
            chunk = file.read(chunk_size)
            for line in chunk.split("\n"):
                values = line.strip().split("\t")
                molecule = dict(zip(keys, values))
                batch.append(molecule)
                if len(batch) >= max_batch_size:
                    working_db.insert(COLLECTION, batch)
                    batch = []

    finished_files_path = directory + 'finished_files'
    with open(finished_files_path, 'a') as f:
        f.write(filename + '\n')
    print(f'Finished moving {filename} to local database')
    if remove:
        os.remove(filename)


def fetch_pdbbind_old(email='samuelmar@gmail.com', pw='J%407iD6i7Lhnugqj'):
    """"""

    # first we log in
    login_url = "http://www.pdbbind.org.cn/ajaxcode.php?action=login"
    payload = {
        "email": email,
        "passwd": pw,
    }
    session = requests.Session()
    response = session.post(login_url, data=payload)
    if response.ok:
        print("Login to pdbbind successful.")
    else:
        print("Login failed. Please check your credentials or the request URL.")

    # now the downloads to local memory
    download_urls = [
        "https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v2020_other_PL.tar.gz",
        "https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v2020_refined.tar.gz"
    ]

    # specify the directory for saving the downloaded files
    script_dir = os.path.dirname(os.path.realpath(__file__))
    save_dir = os.path.join(script_dir, 'pdbbind')
    os.makedirs(save_dir, exist_ok=True)  # Create directory if it doesn't exist

    # download files
    for url in download_urls:
        file_name = url.split("/")[-1]
        save_path = os.path.join(save_dir, file_name)
        
        print(f"Downloading {file_name}...")
        with session.get(url, stream=True) as response:
            response.raise_for_status()
            with open(save_path, 'wb') as out_file:
                for chunk in response.iter_content(chunk_size=8192):
                    out_file.write(chunk)
        print(f"Downloaded {file_name}.")

        print(f"Extracting {file_name}...")
        with tarfile.open(save_path, 'r:gz') as tar:
            tar.extractall(path=save_dir)
        print(f"Extracted {file_name}.")

def fetch_pdbbind(email='samuelmar@gmail.com', pw='J@7iD6i7Lhnugqj', db='universe', collection_name='pdbbind'):
    import time
    import csv
    db = Database(db)
    # first we log in
    login_url = "http://www.pdbbind.org.cn/ajaxcode.php?action=login"
    payload = {
        "email": email,
        "passwd": pw,
    }
    headers = {
        'User-Agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:15.0) Gecko/20100101 Firefox/15.0.1',
        'Referer': 'http://www.pdbbind.org.cn/index.php',
        'X-Requested-With': 'XMLHttpRequest',
    }

    session = requests.Session()
    session.headers.update(headers)
    response = session.post(login_url, data=payload)
    #print(f"Login response: {response.text}")
    response_content = response.json()
    if response_content.get('status') == 'success':
        print("Login to pdbbind successful.")
    else:
        print(f"Login failed. Server response: {response_content.get('message')}")
        return


    # now the downloads to local memory
    download_urls = [
        "https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v2020_other_PL.tar.gz",
        "https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v2020_refined.tar.gz"
    ]

    # specify the directory for saving the downloaded files
    script_dir = os.path.dirname(os.path.realpath(__file__))
    save_dir = os.path.join(script_dir, 'pdbbind')
    os.makedirs(save_dir, exist_ok=True)  # Create directory if it doesn't exist

    block_size = 8192
    # download files
    for url in download_urls:
        file_name = url.split("/")[-1]
        save_path = os.path.join(save_dir, file_name)

        # Check if the file already exists and get the already downloaded size
        if os.path.exists(save_path):
            initial_pos = os.path.getsize(save_path)
        else:
            initial_pos = 0
        
       # Retry download if connection breaks
        retry_counter = 0
        while retry_counter < 100:
            try:
                print(f"Downloading {file_name}...")

                # Resume download if the file already exists
                resume_header = {'Range': f'bytes={initial_pos}-'} if initial_pos else None
                response = session.get(url, headers=resume_header, stream=True)
                if response.status_code not in [200, 206]:  # 200: OK, 206: Partial Content
                    print(f"Download response: {response.status_code}")
                    response.raise_for_status()

                # Get the total size
                if response.status_code == 206:  # Partial Content
                    content_range = response.headers.get('content-range')
                    total_size = int(content_range.split('/')[-1])  # Get the total size from the Content-Range header
                else:
                    total_size = int(response.headers.get('content-length', 0))

                # If the existing file's size is equal to or larger than total_size, skip the download
                if initial_pos >= total_size:
                    print(f"{file_name} is already downloaded.")
                    break

                t = tqdm(total=total_size, initial=initial_pos, unit='iB', unit_scale=True)

                # Write to the file
                with open(save_path, 'ab' if initial_pos else 'wb') as out_file:
                    for data in response.iter_content(block_size):
                        t.update(len(data))
                        out_file.write(data)
                t.close()

                # Verifying downloaded size
                if os.path.getsize(save_path) != total_size:
                    print("ERROR, something went wrong in the downloading process")

                print(f"Downloaded {file_name}.")
                break  # exit the while loop

            except (requests.exceptions.RequestException, requests.exceptions.ChunkedEncodingError) as e:
                print(f"Error occurred while downloading {file_name}: {str(e)}")
                retry_counter += 1
                print(f"Retry attempt: {retry_counter}")
                if os.path.exists(save_path):
                    initial_pos = os.path.getsize(save_path)
                time.sleep(5)  # backoff time before retrying

    ############ AFTER DOWNLOADING, PARSE THEN SEND TO LOCAL DB ############

    extract = False # line for debugging purposes
    if extract:
        # Loop over all tar files in the directory
        for tar_filename in os.listdir(save_dir):
            print(f'Decompressing file: {tar_filename}, this can take a while...')
            if tar_filename.endswith('.tar.gz'):
                # Extract the tar file
                tar_filepath = os.path.join(save_dir, tar_filename)
                tar = tarfile.open(tar_filepath)
                tar.extractall(path=save_dir)
                tar.close()

    insert = True # another line for debugging
    if insert:
        excluded_folders = ['readme']
        file_types = ['.mol2', '.sdf', '.pdb', '.2020']
        used_collections = []
        print('Organizing extracted data, this can take a while...')
        for item in os.listdir(save_dir):
            item_path = os.path.join(save_dir, item)
            # check if the item is a directory
            if os.path.isdir(item_path):
                folder_path = item_path
                print(f'Organizing: {folder_path}')
                collection_name = item
                used_collections.append(collection_name) # save collection names for .2020 files data extracting later
                documents = [] # new list of docs for each collection
                indexinfo = [] # empty list to hold contents of index files
                # Loop over all subfolders (each representing an index)
                for index_name in tqdm(os.listdir(folder_path)):
                    if index_name not in excluded_folders:
                        index_path = os.path.join(folder_path, index_name)

                        if index_name != 'index':
                            # Prepare a document for this index
                            document = {"index": index_name}

                        # Loop over all files in the subfolder
                        for filename in os.listdir(index_path):
                            file_path = os.path.join(index_path, filename)
                            file_ext = os.path.splitext(filename)[1]

                            # If it's not a file or not the file type we are interested in, skip
                            if not os.path.isfile(file_path) or file_ext not in file_types:
                                continue
                            
                            with open(file_path, 'r') as f:
                                file_content = f.read()
  
                            # Add file content to the document under appropriate key
                            if file_ext == '.mol2':
                                document['mol2'] = file_content
                            elif file_ext == '.sdf':
                                document['sdf'] = file_content
                            elif '_pocket.pdb' in filename:
                                document['pocket pdb'] = file_content
                            elif '_protein.pdb' in filename:
                                document['protein pdb'] = file_content
                            elif file_ext == '.2020':
                                indexinfo.append(file_content)
                        if index_name != 'index':
                            documents.append(document)

                # Insert all documents from single tar.gz file into relative MongoDB collection
                print(f'Inserting data into the {collection_name} mongo collection...')
                db.insert(collection_name, documents)

        print(f'Collecting additional data from .2020 index files...')
        # now for the extra data in the index files
        data = {}
        for file2020 in indexinfo:
            # grab only column names from header
            header_line = re.search(r'# PDB code(.*?)\n', file2020).group(1)
            header_line = 'PDB code' + header_line
            column_names = [word.strip() for word in header_line.split(',')]
            column_names.append('extra info') # gotta add this manually because pdbbind is so well made by super smart people

            # Check if we have either 'binding data' or 'Kd/Ki' in column_names
            target_column_names = ['binding data', 'Kd/Ki']
            target_column_index = None
            for target in target_column_names:
                if target in column_names:
                    target_column_index = column_names.index(target)
                    break

            # grab everything but header
            all_lines = file2020.split('\n')
            data_lines = [line for line in all_lines if not line.startswith('#') and line.strip() != '']
            for line in data_lines:
                # each row must be manually parsed because the people at pdbbind are great and intelligent
                line = line.replace('//', '').strip()
                parts = line.split()

                # Split the 'binding data' or 'Kd/Ki' into 'standard type', 'standard relation' and 'standard value' only if the column exists
                if target_column_index is not None:
                    if '=' in parts[target_column_index]:
                        standard_type, standard_value = parts[target_column_index].split('=')
                        standard_relation = '='
                    elif '>' in parts[target_column_index]:
                        standard_type, standard_value = parts[target_column_index].split('>')
                        standard_relation = '>'
                    elif '<' in parts[target_column_index]:
                        standard_type, standard_value = parts[target_column_index].split('<')
                        standard_relation = '<'
                    else:
                        standard_type = parts[target_column_index]
                        standard_relation = 'NA'
                        standard_value = 'NA'
                    parts[target_column_index] = standard_type
                    parts.insert(target_column_index + 1, standard_relation)
                    parts.insert(target_column_index + 2, standard_value)

                    # Update column_names to accommodate the new columns
                    column_names[target_column_index] = 'standard type'
                    column_names.insert(target_column_index + 1, 'standard relation')
                    column_names.insert(target_column_index + 2, 'standard value')

                # assign column values to column names
                row = {column_names[i]: value for i, value in enumerate(parts[:len(column_names) - 1])}

                # If there is extra info, join it into a single string and add it to the dictionary
                if len(parts) > len(column_names) - 1:
                    row['extra info'] = ' '.join(parts[len(column_names) - 1:])

                # Add the row to the data dictionary, using the PDB code as the key
                pdb_code = row['PDB code']
                data[pdb_code] = row

        # reshape additional .2020 data to fit mongo standards, then send to mongo
        documents = [{"index": k, **v} for k, v in data.items()]

        ghost_indexes = []  # list to hold ghost indexes
        print('Sending data to local mongo db...')
        for document in tqdm(documents):
            index_found = False
            for collection_name in used_collections:
                # check if the index exists in the collection
                if db.exists(collection_name, document['index']):
                    db.insert(collection_name, [document])  # update the document
                    index_found = True
                    break  # no need to check other collections
            if not index_found:
                ghost_indexes.append(document['index'])  # save ghost index

        # print out ghost indexes if any
        #if ghost_indexes:
            #print('Ghost indexes found:', ghost_indexes)
            





    