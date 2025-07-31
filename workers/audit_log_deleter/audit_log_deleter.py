import datetime as dt
import logging
import os


def delete_old_messages_from_db(max_lifetime, connection_str: str):
    pass


def main():
    logger = logging.Logger("audit_log_deleter", level=logging.INFO)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    max_lifetime = os.getenv("MAX_LIFETIME")
    db_connection_str = os.getenv(
        "DB_CONNECTION_STRING",
        "postgresql://postgres:postgres@localhost:5432/ai4gauditlog",
    )

    running = True
    while running:
        try:
            logger.info(f"Deleting messages to DB...")
            delete_old_messages_from_db(max_lifetime, db_connection_str)
        except KeyboardInterrupt:
            running = False
            logger.info("Exitiing...")
        except Exception as e:
            running = False
            logger.error(f"An error occurred...{os.linesep}{e}")


if __name__ == "__main__":
    main()
