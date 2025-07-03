FROM quay.io/minio/minio:RELEASE.2025-04-22T22-12-26Z

# Copy the custom entrypoint script
COPY minio-setup.sh /usr/local/bin/minio-entrypoint.sh
RUN chmod +x /usr/local/bin/minio-entrypoint.sh

# Use the custom entrypoint
ENTRYPOINT ["/usr/local/bin/minio-entrypoint.sh"]

# Default CMD to run MinIO server
CMD ["server"]
