#About
Orrery is a Python-based project designed to perform molecular parsing and data handling with MongoDB. The project includes a Docker setup for easy deployment and isolation of dependencies.

#Basic Usage
1. Ensure Docker and Docker-compose are installed on your machine.

2.Navigate into the orrery directory:
    ~ cd orrery

3.Build and start the Docker containers:
    ~ docker-compose up -d

    or, in the event of rebuilding, first run:
    ~ docker-compose build
    
4. Access a bash terminal inside the desired container:
    ~ docker exec -it TESTer <container_name> bash

5. Execute commands normally as you would in any virtual environment.

#Modules
classes.py: Contains the definition of the Molecule class used throughout the project.
parsers.py: Provides functionality for parsing mol2 files.
mongo_api.py: Provides a range of functions for interacting with MongoDB.

Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

License
MIT