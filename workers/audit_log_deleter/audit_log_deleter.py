import datetime as dt
import logging
import os

from sqlalchemy import Engine, create_engine, text
from sqlalchemy.orm import Session


def delete_old_events_from_db(engine: Engine, n_days: int):
    """Delete old audit log events from the database. Events older than
    today - n_days will be deleted.

    Args:
        engine (Engine): The engine that connects to the database.
        n_days (int): Messages older than this will be deleted.
    """
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

    max_lifetime = int(os.getenv("MAX_LIFETIME", 1))
    db_connection_str = os.getenv(
        "DB_CONNECTION_STRING",
        "postgresql://postgres:postgres@localhost:5434/ai4gauditlog",
    )
    # Establish DB connection
    engine = create_engine(db_connection_str)

    running = True
    while running:
        try:
            delete_old_events_from_db(engine, max_lifetime)
        except KeyboardInterrupt:
            running = False
            logger.info("Exitiing...")
        except Exception as e:
            running = False
            logger.error(f"An error occurred...{os.linesep}{e}")


if __name__ == "__main__":
    main()
