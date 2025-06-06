# AI4Green Installation Guide
___

AI4Green is an Electronic Laboratory Notebook which combines data storage and sharing while promoting green and sustainable chemistry.

**We highly recommend using Chromium-based browsers such as Google Chrome for accessing AI4Green.**

See the [video tutorials](https://youtube.com/playlist?list=PL7u_tOd0vTynC0tgWt0cb7jLlnEoDwiO4) on how to use AI4Green.

___
## Option 1: _Via_ Our Web Application

Visit [ai4green.app](https://ai4green.app/) for more information and to register for a free account.

___

## Option 2: Build from source _via_ Conda

### Step 1 - Download Source Code

The source code can be downloaded from this repository. Place the contents in a folder named `AI4Green`.
Windows users should run from following commands from Anaconda Prompt and Mac OS/Linux users from terminal.

### Step 2 - Set up Virtual Environment

We recommend using Anaconda to set up the python environment, available at [anaconda.com](https://www.anaconda.com/).

`cd AI4Green`

`conda create --name ai4green python==3.10`

`conda activate ai4green`

`pip install poetry==1.8.3`

`poetry install`

This will install all the dependencies listed in pyproject.toml into the current conda environment.

Error with lib magic? [See here](#python-magic)

### Step 3 - App Configuration

Configuration of AI4Green is controlled by `AI4Green/Webapp/sources/config.py`, the default values are included, but can be overridden by setting environment variables.

`SQLALCHEMY_DATABASE_URI`: Postgres database connection string. 

`COMPOUND_LIMIT`: Sets the limit for the number of compounds extracted from PubChem. Note that a large database can cause smaller servers to crash.

`MAIL`: Ai4Green is configured to send email alerts triggered by actions such as managing and creating Workgroups. Functionality such as password reset requires email support. 
See `AI4Green/Webapp/sources/config.py` lines `55 to 64` to edit this configuration.

We use `Flask-Mail`, which can be set up by providing the email address and password. Note that you may need to configure your email account to allow messages to be sent from the application.

`MARVIN_JS_API_KEY`: API key for Marvin JS is specified here. If not provided an alternative sketcher is available.

### Step 4 - Build the Database

The AI4Green database can be built with docker. This required docker to be installed on your system. From the AI4Green directory run

`docker-compose up -d`

Next, with your conda environment activated, navigate to the Webapp directory

`cd Webapp`

Make the database migrations 

`flask db upgrade`

and download the pubchem database

`flask download-pubchem`

The user is prompted to select whether to use an existing download or PubChem database or download the latest version.

Finally, the user decides whether to add PubChem Compound data to the existing database (i.e. retain User, Workgroup, Workbook and Reaction data) or create a new database.
`flask download-pubchem` can be rerun at any time to update the database.

To seed the data from the pubchem download, run

`flask seed-db` 

To update your version of pubchem at any time, run

`flask update-pubchem`

#### Seed data (optional)

Seed the data in `seed_data.yaml` to include default users

`flask seed-users` 

Otherwise, you can create a new user through the CLI, by running 

`flask add-user` 

and following the prompts.

### Step 5 - Run AI4Green

The flask app can now be run locally.

`conda activate ai4green`

`cd AI4Green/Webapp`

`flask run`

AI4Green can now be accessed via your browser at [127.0.0.1:80/home/](http://127.0.0.1:80/home).

### Step 6 - Deploy to Server (Optional)

It is possible to deploy an instance of AI4Green to a remote server for increased data security and privacy. 

A system administrator should be able to set up a version which is available on an organisation's internal server.

There are several useful guides on how to deploy a web application with a Postgres database on cloud services For example, Microsoft Azure provide a collection of [useful guides](https://learn.microsoft.com/en-us/azure/?product=popular).
___
## Development

For developers running the application locally.

### Prerequesites

- Poetry
- Docker

### Getting Started
This is similar to the default setup, but you must also install the pre-commit hooks:
1. Install the dependencies: `poetry install`
2. Install pre-commit hooks: `pre-commit install`
3. Use the CLI to update the database, follow [step 4](#step-4---build-the-database)
4. Run the app `flask run`

### Database Migrations

Database changes are handled using [Flask-Migrate](https://flask-migrate.readthedocs.io/en/latest/#) and you should read their documentation for usage.

Briefly, to create a new migration after you have made changes in `sources/models/`

```bash
flask db upgrade # ensure current database is up to date
flask db migrate -m "Migration name" # create the migration
flask db upgrade # update the database
```

### Apache Kafka
The `docker-compose.yaml` file contains the `kafka` service (which is supported by `zookerper`) and the `kafka-ui` service. `kafka-ui` can be used to set up a cluster and topics in the browser at http://localhost:8080. 

**When using Kafka, make sure you set the environment variable `USE_KAFKA` to 1. Otherwise messages will be logged to the standard output. This behaviour is only for testing purposes where a Kafka instance is not required.**

When you first start the services and navigate to the `kafka-ui` in the browser, you will be asked to configure a cluster. Enter a name in the "Cluster name" field, then under "Bootstrap Servers", enter "kafka" in the "Host" field and 29092 in the "Port" field. Click the "Validate" button to make sure it's working, then click "Submit" if you get the green dialog box.

Next you can start creating topics. On the left-side of the screen, you will see a menu with your configured cluster(s) below "Dashboard". Click your cluster to expand the drop-down and click "Topics". In the top-right, click the "Add a Topic" button. Give your topic a name, choose a number of partitions, the cleanup policy (you can use delete for testing purposes) and choose a retention time (there are presets to choose from). Then press "Create topic".

When running the flask app, you can set the `MESSAGE_QUEUE_HOSTNAME` environment variable in a `.env` file or it defaults to `localhost:9092`.

### Kafka Connect + MinIO for long term log persistence
To store event logs long-term, there are 2 more services in `docker-compose.yaml` called `kafka-connect` and `minio`. `kafka-connect` acts as a consumer of events from `kafka` and can redirect them to other services, in this case `minio`. `minio` is an object storage solution that is compatible with AWS S3 APIs. `kafka-connect` consumes the events generated by AI4Green and stores them long term in `minio` for audit purposes.

First, set up `minio`. Open `docker-compose.yaml` and go to the `minio` service. Set `MINIO_ROOT_USER` and `MINIO_ROOT_PASSWORD` to values of your choice and run the service. With the service running from `docker-compose.yaml`, navigate to http://localhost:9001 and sign in. Create a bucket for your kafka events and note the name. Next create an access key for `minio` and save the ID/Secret pair.

To set up Kafka Connect, you will first need to download the [Amazon S3 Sink Connector](https://www.confluent.io/hub/confluentinc/kafka-connect-s3). Click the "Download" button under "Self-Hosted". You will get a `.zip` file with the plugin. Unzip the plugin and place the unzipped contents into a directory called `connect-plugins`, which you should make in same directory as `docker-compose.yaml`. Now run the `kafka-connect` service.

Once the 2 services are running, you can connect `kafka-connect` to `kafka` using `kafka-ui`. Navigate to http://localhost:8080 to view your `kafka` topics. Click "Configure" on the topic you wish to consume with `kafka-connect`. Scroll down to where it says "Kafka Connect" and click the button saying "Configure Kafka Connect". Give the connection a name and put "http://kafka-connect:8083" in the field marked "Kafka Connect URL". Then, click "Validate" at the buttom of the page to make sure the information is valid, and then press "Submit".

Now that `kafka-connect` is linked to your topic in `kafka`, you need to configure a connector to consume the events from your topic. In the directory `kafka-connect-configs` is a file called `example.json`. You can copy this file and replace various fields to configure your connector.
Replace the value in `"topics"` with the topic(s) you wish to consume, either as a single value or a comma-separated string of values. Replace `"s3.bucket.name"` with the name you gave your `minio` bucket, `"aws.access.key.id"` with the ID of the access key you created in `minio`, and `"aws.secret.access.key"` with the Secret. You should also set `"flush.size"` to a larger number, especially in production. This number determines how many messages need to be sent to your topic before `kafka-connect` writes to `minio`. A lower number is useful for test purposes, though. Update `"partition.field.name"` to either a single value or comma-separated string of values. The order of values determines the structure of the folders in `minio`.

Once you have finished your configuration `.json` file, you need to tell `kafka-connect` to register it. You can use `curl` in the terminal. From the `kafka-connect-configs` directory, run:

```bash
curl -X POST -H "Content-Type: application/json" --data @<your config>.json http://localhost:8083/connectors
```
___
## Troubleshooting

### Conda SSL Module Error

- An error with Anaconda on Windows 10 has been reported, when using the `conda install` command. Details of this error and a solution can be found [here](https://github.com/conda/conda/issues/8273).


- An error has been reported using `poetry install` on Windows 10 and 11 machines. Using `poetry==1.8.3` should fix this issue.

### Error Creating Database

- It is a known issue that PubChem sometimes change the format of their data between the releases, which causes an error extracting that data.
A previous version of the PubChem compound data can be obtained by emailing [ai4green@nottingham.ac.uk](mailto:ai4green@nottingham.ac.uk).

### Marvin JS Sketcher Error

- If there is a problem with Marvin JS or you are deploying AI4Green to another url rather than running it locally, you may need a new Marvin JS API key and register the domain. You can find these at [pro.chemicalize.com/app](https://pro.chemicalize.com/app/marvin/settings). Please note you will have to register for a ChemAxon account to use this service. The Marvin JS api key can then be replaced in `AI4Green/Webapp/config.py`.


### Python Magic

- There have been reports of python magic not being installed correctly using poetry, giving this error

`ImportError: failed to find libmagic. Check your installation`

- This can often be fixed by

`pip install python-magic-bin`

