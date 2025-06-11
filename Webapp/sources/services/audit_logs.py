from minio import Minio
from minio.error import S3Error


def main():
    # Create a client with the MinIO server playground, its access key
    # and secret key.
    client = Minio(
        "localhost:9000",
        access_key="wwBrRE5ztVSHaZiFpFN0",
        secret_key="Lfbw26l1uoiR6FQ4JFZPOMWWI8LJ7lLbQnQuuzIO",
        secure=False,
    )

    # The destination bucket and filename on the MinIO server
    bucket_name = "kafka-logs"

    # Make the bucket if it doesn't exist.
    found = client.bucket_exists(bucket_name)
    if not found:
        client.make_bucket(bucket_name)
        print("Created bucket", bucket_name)
    else:
        print("Bucket", bucket_name, "already exists")

    # Upload the file, renaming it in the process
    print(client.list_buckets())


if __name__ == "__main__":
    try:
        main()
    except S3Error as exc:
        print("error occurred.", exc)
