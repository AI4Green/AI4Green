from enum import Enum

# all enums for ai4green app are stored here


class ReactionType(Enum):
    """
    Enum for labelling reaction types, used on reaction save
    """

    STANDARD = "STANDARD"
    POLYMER = "POLYMER"


class ApprovalStatus(Enum):
    """
    Enum for labelling approval statuses for data export requests
    """

    PENDING = "PENDING"
    APPROVED = "APPROVED"
    REJECTED = "REJECTED"
    EXPIRED = "EXPIRED"


class ExportFormat(Enum):
    """
    Enum for labelling export formats for data export requests
    """

    RDF = "RDF"
    RXN = "RXN"
    PDF = "PDF"
    ELN = "ELN"
    SURF = "SURF"
    CSV = "CSV"
    JSON = "JSON"
    SI = "SI"


class ReactionApprovalStatus(Enum):
    """
    Enum for labelling reaction approval statuses for reaction review process
    """

    PENDING = "PENDING"
    APPROVED = "APPROVED"
    REJECTED = "REJECTED"
    CHANGES_REQUESTED = "CHANGES_REQUESTED"


class SubscriptionTier(Enum):
    """
    Enum for tracking subscription tiers for workgroups
    """

    FREE = "FREE"
    PRO = "PRO"
    ENTERPRISE = "ENTERPRISE"
