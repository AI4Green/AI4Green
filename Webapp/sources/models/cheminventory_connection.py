from sources.extensions import db


class ChemInventoryConnection(db.Model):
    """A workgroup's link to its ChemInventory account. One per workgroup.

    The API token is stored encrypted (services.cheminventory.encrypt_token) and
    is never sent back to the browser. `token_hint` holds its last few
    characters so people can tell which token is in use.
    """

    __tablename__ = "ChemInventoryConnection"

    id = db.Column(db.Integer, primary_key=True)
    workgroup_id = db.Column(
        db.Integer, db.ForeignKey("WorkGroup.id"), nullable=False, unique=True
    )
    encrypted_token = db.Column(db.Text, nullable=False)
    token_hint = db.Column(db.String(8), nullable=False)

    # From ChemInventory when the token was checked, for display only
    cheminventory_user = db.Column(db.String(255))
    inventory_name = db.Column(db.String(255))

    created_by_id = db.Column(db.Integer, db.ForeignKey("User.id"), nullable=False)
    created_at = db.Column(
        db.DateTime(timezone=True), default=db.func.now(), nullable=False
    )
    updated_at = db.Column(
        db.DateTime(timezone=True),
        default=db.func.now(),
        onupdate=db.func.now(),
        nullable=False,
    )

    workgroup = db.relationship(
        "WorkGroup",
        backref=db.backref(
            "cheminventory_connection", uselist=False, cascade="all, delete-orphan"
        ),
    )
    created_by = db.relationship("User")
