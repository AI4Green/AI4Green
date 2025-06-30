#!/bin/bash

set -e

# Wait for Kafka Connect to be available
echo "Waiting for Kafka Connect to be available at $CONNECT_REST_URL..."
until curl -s "$CONNECT_REST_URL/connectors" > /dev/null; do
  sleep 3
done

echo "Kafka Connect is available. Creating connector..."

# Build JSON config dynamically using env variables
cat <<EOF > /tmp/connector-config.json
{
  "name": "minio-sink-connector",
  "config": {
    "topics": "${KAFKA_TOPICS}",
    "s3.bucket.name": "${S3_BUCKET_NAME}",
    "partition.field.name": "${PARTION_FIELDS}",
    "store.url": "${MINIO_URL}",
    "s3.endpoint": "${MINIO_URL}",
    "aws.access.key.id": "${AWS_ACCESS_KEY_ID}",
    "aws.secret.access.key": "${AWS_SECRET_ACCESS_KEY}",
    "connector.class": "io.confluent.connect.s3.S3SinkConnector",
    "tasks.max": "1",
    "s3.region": "us-east-1",
    "s3.part.size": "5242880",
    "flush.size": "3",
    "storage.class": "io.confluent.connect.s3.storage.S3Storage",
    "format.class": "io.confluent.connect.s3.format.json.JsonFormat",
    "schema.compatibility": "NONE",
    "s3.path.style.access": "true",
    "partitioner.class": "io.confluent.connect.storage.partitioner.FieldPartitioner",
    "value.converter": "org.apache.kafka.connect.json.JsonConverter",
    "value.converter.schemas.enable": "true",
    "key.converter": "org.apache.kafka.connect.storage.StringConverter"
  }
}
EOF

# POST the connector configuration
curl -X POST -H "Content-Type: application/json" --data @/tmp/connector-config.json "$CONNECT_REST_URL/connectors"

# Start Kafka Connect normally
echo "Starting Kafka Connect..."
exec /etc/confluent/docker/run
