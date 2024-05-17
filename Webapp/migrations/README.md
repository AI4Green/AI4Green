## Single-database configuration for Flask.

Generated migration scripts may need to be edited to enable commands flask db downgrade and flask db upgrade to work.
This was found using alembic 1.13.1 and may change in future versions.

Manual edits should follow the naming convention defined in the documentation:

```
convention = {
    "ix": 'ix_%(column_0_label)s',
    "uq": "uq_%(table_name)s_%(column_0_name)s",
    "ck": "ck_%(table_name)s_%(constraint_name)s",
    "fk": "fk_%(table_name)s_%(column_0_name)s_%(referred_table_name)s",
    "pk": "pk_%(table_name)s"
}
```

https://flask-sqlalchemy.palletsprojects.com/en/2.x/config/#using-custom-metadata-and-naming-conventions


#### Example migration: 57cdfc52c9ed_data_export_request
This migration script has been edited to add a name to the foreign key constraint.\
The foreign key constraint cannot be dropped if it does not have a name, giving an error during flask db downgrade. \
\
Original:
```
    batch_op.create_foreign_key(
        None,
        "DataExportRequest",
        ["data_export_request_id"],
        ["id"],
    )
```

After-editing:
```
    batch_op.create_foreign_key(
        "fk_Reactions_data_export_request_ids_DataExportRequests",
        "DataExportRequest",
        ["data_export_request_id"],
        ["id"],
    )
```
The downgrade function has been edited to match this by naming the foreign key constraint being dropped and
to include dropping of the ENUMs which are added during the upgrade. \
\

Original:
```
    batch_op.drop_constraint(
        None,
        type_="foreignkey",
    )
```
After-editing:
```
    batch_op.drop_constraint(
        "fk_Reactions_data_export_request_ids_DataExportRequests",
        type_="foreignkey",
    )

    op.execute("DROP TYPE IF EXISTS approvalstatus;")
    op.execute("DROP TYPE IF EXISTS exportformat;")
```
