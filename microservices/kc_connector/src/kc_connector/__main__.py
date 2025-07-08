import json
import logging
import os
import sys
import requests

KAFKA_CONNECT_URL = os.getenv("KAFKA_CONNECT_URL") or exit(
    "KAFKA_CONNECT_URL is not set"
)
MINIO_URL = os.getenv("MINIO_URL") or exit("MINIO_URL is not set")
MINIO_BUCKET_NAME = os.getenv("MINIO_BUCKET_NAME") or exit(
    "MINIO_BUCKET_NAME is not set"
)
MINIO_ACCESS_KEY = os.getenv("MINIO_ACCESS_KEY") or exit("MINIO_ACCESS_KEY is not set")
MINIO_SECRET_KEY = os.getenv("MINIO_SECRET_KEY") or exit("MINIO_SECRET_KEY is not set")
TOPIC = os.getenv("TOPIC") or exit("TOPIC is not set")
PARTITION_FIELDS = os.getenv("PARTITION_FIELDS") or exit("PARTITION_FIELDS is not set")
FLUSH_SIZE = os.getenv("FLUSH_SIZE") or exit("FLUSH_SIZE is not set")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # get the connector name from the command line
    connector_name = sys.argv[1]

    # load the template JSON
    with open("template.json", "r") as template_file:
        template = json.load(template_file)

    # Fill in the configuration
    template["name"] = connector_name
    template["config"]["topics"] = TOPIC
    template["config"]["s3.bucket.name"] = MINIO_BUCKET_NAME
    template["config"]["partition.field.name"] = PARTITION_FIELDS
    template["config"]["store.url"] = MINIO_URL
    template["config"]["s3.endpoint"] = MINIO_URL
    template["config"]["aws.access.key.id"] = MINIO_ACCESS_KEY
    template["config"]["aws.secret.access.key"] = MINIO_SECRET_KEY
    template["config"]["flush.size"] = FLUSH_SIZE

    # Send the data to kafka-connect
    res = requests.post(f"{KAFKA_CONNECT_URL}/connectors", json=template)
    try:
        res.raise_for_status()
        logging.info(f"Successfully created connector {connector_name}")
    except:
        logging.error(
            f"Failed to create connector {connector_name}: {res.status_code} {res.reason}"
        )
    finally:
        res.close()

    # Clean up the logger
    logging.shutdown()
