import json
from typing import Any, Callable, List
from minio import Minio
from minio.error import S3Error


def _extract_logs(logs: list, deserialiser: Callable) -> List[Any]:
    """Decode the logs from a S3 storage object.

    Args:
        logs (list): The raw list of logs from the S3 API containing the log records.
        deserialiser (Callable):
            The function to deserialise the log records. This should be the `deserialise`
            static method on a message class with the `MessageSerdeMixin`.

    Returns:
        List[Any]: The log records.
    """
    extracted = []
    for line in logs:
        extracted.append(deserialiser(json.loads(line)))
    return extracted


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
    objects = list(client.list_objects(bucket_name, recursive=True))
    response = client.get_object(bucket_name, objects[0].object_name)
    if response:
        # print(
        #     *_extract_logs(response.readlines(), ReactionEditMessage.deserialise),
        #     sep="\n"
        # )
        response.close()


if __name__ == "__main__":
    try:
        main()
    except S3Error as exc:
        print("error occurred.", exc)
