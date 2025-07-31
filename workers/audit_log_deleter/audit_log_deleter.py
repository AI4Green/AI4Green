import datetime as dt
import logging
import os
from sqlalchemy import create_engine, text
from sqlalchemy.orm import Session


def delete_old_events_from_db(n_days: int, connection_str: str):
    """Delete old audit log events from the database. Events older than
    today - n_days will be deleted.

    Args:
        n_days (int): Messages older than this will be deleted.
        connection_str (str): The connection string to the database.
    """
    # Establish DB connection
    engine = create_engine(connection_str)
    time_limit = dt.datetime.now() - dt.timedelta(days=n_days)

    # Create the session
    with Session(engine) as session:
        session.begin()
        try:
            # delete old records
            query = text('DELETE FROM "AuditLogEvent" WHERE event_time < :time_limit')
            session.execute(query, {"time_limit": time_limit})
        except:
            # rollback if anything goes wrong
            session.rollback()
            raise
        else:
            # if all goes well, commit the changes
            session.commit()


def main():
    logger = logging.Logger("audit_log_deleter", level=logging.INFO)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    max_lifetime = os.getenv("MAX_LIFETIME", 1)
    db_connection_str = os.getenv(
        "DB_CONNECTION_STRING",
        "postgresql://postgres:postgres@localhost:5432/ai4gauditlog",
    )

    running = True
    while running:
        try:
            logger.info(f"Deleting messages to DB...")
            delete_old_events_from_db(max_lifetime, db_connection_str)
        except KeyboardInterrupt:
            running = False
            logger.info("Exitiing...")
        except Exception as e:
            running = False
            logger.error(f"An error occurred...{os.linesep}{e}")


if __name__ == "__main__":
    main()
