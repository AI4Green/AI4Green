from flask_login import current_user
from typing import List
from sources import models
from sources.extensions import db


def list_user_graphs(user: models.Person) -> List[models.PCAGraph]:
    """
        Retrieve a list of PCA graphs associated with the current user ordered by date

        Returns:
            List of saved PCA_graphs with most recent first
    """
    graph_query = (
        db.session.query(models.PCAGraph)
        .filter(models.PCAGraph.creator_person == user)
        .filter(models.PCAGraph.status == "active")
    )

    return (
        graph_query.order_by(models.PCAGraph.time_of_creation.desc())
        .all()
    )


def PCAgraph_from_id(graph_id: str, user: models.Person) -> models.PCAGraph:

    """
        Retrieve the PCA_graph corresponding to graph_id associated with the current user

        Returns:
            models.PCAGraph with corresponding graph_id
    """

    return (
        db.session.query(models.PCAGraph)
        .filter(models.PCAGraph.creator_person == user)
        .filter(models.PCAGraph.id == graph_id)
        .first()
    )


def list_all() -> List[models.PCAGraph]:
    """
        Gets a list of all PCA graphs saved. For the admin_dashboard

        Returns:
             List of all PCA graphs sorted by date
        """
    return (db.session.query(models.PCAGraph)
            .order_by(models.PCAGraph.time_of_creation.desc())
            ).all()
