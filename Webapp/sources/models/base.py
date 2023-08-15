# coding: utf-8
from sources.extensions import db
from typing import Any

class CRUDMixin(object):
    """
    Mixin that adds convenience methods for CRUD (create, read, update, delete) operations.
    """

    @classmethod
    def create(cls, commit: bool = True, **kwargs: Any) -> object:
        """Create a new record and save it the database."""
        instance = cls(**kwargs)
        return instance.save() if commit else instance.save(commit=False)

    def update(self, commit: bool = True, **kwargs: Any) -> object:
        """Update specific fields of a record."""
        for attr, value in kwargs.items():
            setattr(self, attr, value)
        return self.save() if commit else self

    def save(self, commit: bool = True) -> object:
        """Save the record."""
        db.session.add(self)
        if commit:
            db.session.commit()
        return self

    def delete(self, commit: bool = True) -> None:
        """Remove the record from the database."""
        db.session.delete(self)
        if commit:
            return db.session.commit()
        return


class Model(CRUDMixin, db.Model):
    """
    Base model class that includes CRUD convenience methods.
    """

    __abstract__ = True
