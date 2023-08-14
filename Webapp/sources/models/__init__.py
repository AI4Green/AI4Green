"""

This module contains the models for the database.

The models are used to create the database tables, and to query the database.

"""

from . import base
from .compound_data_error_report import CompoundDataErrorReport
from .compound import Compound
from .element import Element
from .hazard_code import HazardCode
from .institution import Institution
from .news_item import NewsItem
from .notification import Notification
from .novel_compound import NovelCompound
from .person import (
    Person,
    t_Person_WorkBook,
    t_Person_WorkGroup,
    t_Person_WorkGroup_2,
    t_Person_WorkGroup_3,
)
from .reaction_data_file import ReactionDataFile
from .reaction_note import ReactionNote
from .reaction import Reaction, t_Reaction_Reaction
from .role import Role
from .solvent import Solvent
from .user import User
from .wb_status_request import WBStatusRequest
from .wg_status_request import WGStatusRequest
from .workbook import WorkBook
from .workgroup_request import WorkGroupRequest
from .workgroup import WorkGroup
from .update_compound import UpdateCompound
