#!/bin/sh

apt update -y
apt install libxext-dev -y
apt install libxrender1 -y
apt install libmagic-dev -y
gunicorn "Webapp.sources:create_app('prod')" --bind=0.0.0.0 --timeout 600 --chdir ./Webapp
