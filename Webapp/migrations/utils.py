def perform_migration(db_name: str, target: str = "default") -> bool:
    """Create a boolean to determin if a migration should be performed on the
    named database (`db_name`). If `db_name` and `target` match, the migration
    can be performed.

    Args:
        db_name (str): The name of the database to migrate.
        target (str, optional): The intended database to migrate. Defaults to "default".

    Returns:
        bool: Whether the migration can be performed.
    """
    return db_name == target
