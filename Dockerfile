# init a base image
FROM continuumio/miniconda3
# Create the environment:
COPY environment.yaml .
RUN conda env create -f environment.yaml
# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "ai4green", "/bin/bash", "-c"]
# install RDKit
RUN conda install -c conda-forge rdkit
# run pip to install the dependencies of the flask app
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt
# The code to run when container is started:
COPY . .
EXPOSE 80
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "ai4green", "python", "Webapp/app.py"]
