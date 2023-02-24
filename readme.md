# AI4Green Installation Guide
AI4Green is an Electronic Laboratory Notebook which combines data storage and sharing while promoting 
green and sustainable chemistry.<br><br>
<b>We highly recommend using Chromium-based browsers such as Google Chrome for accessing AI4Green.</b><br><br>
See the following video tutorials on how to use AI4Green:
https://youtube.com/playlist?list=PL7u_tOd0vTynC0tgWt0cb7jLlnEoDwiO4


## Option 1: <i>Via</i> Our Web Application
Visit https://ai4green.app/ for more information and to register for a free account.

## Option 2: <i>Via</i> Docker

### Step 1 - Install Docker
Download Docker Desktop for your distribution:<br>
https://docs.docker.com/get-docker/

### Step 2 - Run Docker Image
From your terminal pull and run the ai4green image:<br><br>
`docker run --rm -p 80:80 ai4green/ai4green`<br><br>
An instance of AI4Green is now running locally on your machine.<br>
This can be accessed via your browser at [127.0.0.1:80/home/](). If port 80 is already in use, 
you may change this by running e.g. `-p 2000:80` from Docker then accessing the app at the corresponding port.
<br><br>
Within the docker image there is an admin account with the following credentials:<br><br>
Username: Admin<br>
Password: Admin_login<br><br>
The database is built with SQLite and contains information from 100,000 compounds from PubChem. 
It is not possible to interrogate or update the database using this docker image. Note that this option is primarily 
for evaluation. All data is lost when the container is stopped.

## Option 3: Build from source <i>via</i> Conda

### Step 1 - Download Source Code
The source code can be downloaded from this repository. Place the contents in a folder named `AI4Green`. 
Windows users should run from following commands from Anaconda Prompt and Mac OS/Linux users from terminal.<br>

### Step 2 - Set up Virtual Environment
We recommend using Anaconda to set up the python environment, available at https://www.anaconda.com/. <br><br>
`cd AI4Green`<br>
`conda env create -f environment.yaml`<br>
`conda activate ai4green`<br>
`conda install -c conda-forge rdkit`<br>
`python -m pip install -r requirements.txt`<br>

### Step 3 - Examine Configuration
Configuration of AI4Green is controlled by `AI4Green/Webapp/configs.yaml`. Default settings may be used in the first 
instance, or customised by modifying the following sections. Please note that you must define admin credentials before 
setting up the database.<br><br>
`database_configurations:`<br>
This sets up the database, predefined users, 
and email server for the web application. The database can be configured to use SQLite or Postgres, 
we highly recommend using SQLite in the first instance.<br><br>
`compound_limit:`<br>
This sets the limit for the number of compounds extracted from PubChem. Note that a large database can cause smaller 
servers to crash.<br><br>
`predefined_users:`<br>
You must predefine an admin account before setting up the database. Username, email, password, and fullname 
should be provided.<br><br>
`mail:`<br>
Ai4Green is configured to send email alerts triggered by actions such as managing and creating Workgroups.
Functionality such as password reset requires email support.
We use `flask_mail` which can be set up by providing the email address and password. Note that you may need to configure 
your email account to allow messages to be sent from the application.<br><br>
`marvin_js_api_key:`<br>
API key for Marvin JS is specified here.

### Step 4 - Set up Postgres (Optional)
Information is provided at https://www.postgresql.org/download/ on how to set up Postgres on your distribution.
Note that Window users may have to add Postgres to their system variables. To use PGAdmin, the port, username, 
and password should match the setting in `AI4Green/configs.yaml`.

### Step 5 - Build the Database
Next we run `AI4Green/Compound_database/pubchem_download.py` to create the database. 
The database may be built using SQLite or Postgres. We highly recommend using SQLite initially.<br><br>
`conda activate ai4green`<br>
`cd AI4Green/Webapp/compound_database`<br>
`python pubchem_download.py`<br>

The user is prompted to select whether to use an existing download or PubChem database or download the latest version.
Next the type of database is specified: SQLite or Postgres. Finally, the user decides whether to add PubChem Compound 
data to the existing database (i.e. retain User, Workgroup, Workbook and Reaction data) or create a new database.
`pubchem_download.py` can be rerun at any time to update the database.

### Step 6 - Run AI4Green
The flask app can now be run locally.<br><br>
`conda activate ai4green`<br>
`cd AI4Green/Webapp`<br>
`flask run`<br>

AI4Green can now be accessed via your browser at [127.0.0.1/home/]().

### Step 7 - Deploy to Server (Optional)
It is possible to deploy an instance of AI4Green to a remote server for increased data security and privacy. 
A system administrator should be able to set up a version which is available on an organisation's internal server.
<br><br>
There are several useful guides on how to deploy a web application with a Postgres database on cloud services. 
For example, Microsoft Azure provide a collection of useful guides:<br><br>
https://learn.microsoft.com/en-us/azure/?product=popular

### Troubleshooting
#### Conda SSL Module Error
An error with Anaconda on Windows 10 has been reported, when using the `conda install` command. Details of this error 
and a solution can be found here: 
https://github.com/conda/conda/issues/8273.

#### Error Creating Database
It is a known issue that PubChem sometimes change the format of their data between the releases, which causes an error 
extracting that data. 
A previous version of the PubChem compound data can be obtained by emailing admin@ai4green.app.

#### Marvin JS Sketcher Error
If there is a problem with Marvin JS or you are deploying AI4Green to another url rather than running it locally, 
you may need a new Marvin JS API key and register the domain. You can find these at 
https://pro.chemicalize.com/app/marvin/settings. Please note you will have to register for a ChemAxon account to use 
this service. The Marvin JS api key can then be replaced in `AI4Green/Webapp/configs.yaml`.
