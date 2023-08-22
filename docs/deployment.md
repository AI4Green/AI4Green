# Deployment

When deploying the Flask application Azure App service, your custom startup command needs to point to the `startup.sh` file in the root of the repository.
This script will install the necessary binaries, run database migrations, and run the gunicorn server that will run the application.

[Read about this pattern](
https://learn.microsoft.com/en-us/azure/developer/python/configure-python-web-app-on-app-service)
