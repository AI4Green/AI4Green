# Kafka Connector
This mini-project is used to connect the Kafka Connect instance to Kafka so that messages needed for audit logs can be captured by various topics, and stored where they can be read when needed.

## Usage
The script in `src/kc_connector/__main__.py` is run using the workflow `.github/workflows/make-kc-connectors.yml` either on merges to `main` or on `workflow_dispatch` (manually run in GitHub).

In the workflow, the necessary environment variables are specified. The following should come from github secrets:
- `KAFKA_CONNECT_URL`
- `MINIO_URL`
- `MINIO_BUCKET_NAME`
- `MINIO_ACCESS_KEY`
- `MINIO_SECRET_KEY`

`FLUSH_SIZE` can be modified in the workflow as and when needed.