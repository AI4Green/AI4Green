#!/bin/sh
set -e

# Validate required environment variables
: "${MINIO_ROOT_USER:?Environment variable MINIO_ROOT_USER must be set}"
: "${MINIO_ROOT_PASSWORD:?Environment variable MINIO_ROOT_PASSWORD must be set}"
: "${MINIO_BUCKET_NAME:?Environment variable MINIO_BUCKET_NAME must be set}"
: "${MINIO_CLIENT_ACCESS_KEY:?Environment variable MINIO_CLIENT_ACCESS_KEY must be set}"
: "${MINIO_CLIENT_SECRET_KEY:?Environment variable MINIO_CLIENT_SECRET_KEY must be set}"

# Export MinIO credentials for the server
export MINIO_ROOT_USER
export MINIO_ROOT_PASSWORD

# Start MinIO server in background to configure it
/minio server /data &

# Wait for MinIO to become healthy
echo "Waiting for MinIO to be ready..."
until curl -s http://localhost:9000/minio/health/ready >/dev/null; do
  sleep 1
done

# Configure mc client
mc alias set local http://localhost:9000 "$MINIO_ROOT_USER" "$MINIO_ROOT_PASSWORD"

# Create the bucket if it doesn't exist
if ! mc ls local/"$MINIO_BUCKET_NAME" >/dev/null 2>&1; then
  mc mb local/"$MINIO_BUCKET_NAME"
  echo "Bucket $MINIO_BUCKET_NAME created."
else
  echo "Bucket $MINIO_BUCKET_NAME already exists."
fi

# Create service account
if ! mc admin user svcacct info local "$MINIO_CLIENT_ACCESS_KEY" >/dev/null 2>&1; then
  mc admin user svcacct add local \
    --access-key "$MINIO_CLIENT_ACCESS_KEY" \
    --secret-key "$MINIO_CLIENT_SECRET_KEY" \
    "$MINIO_ROOT_USER"
  echo "Service account $MINIO_CLIENT_ACCESS_KEY created."
else
  echo "Service account $MINIO_CLIENT_ACCESS_KEY already exists."
fi

# Kill background MinIO process
kill %1
wait %1 2>/dev/null || true

# Start MinIO server in foreground
exec /minio server /data
