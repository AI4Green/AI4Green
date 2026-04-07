/**
 * Contains constants for the UI of the experiment.
 */

import {
  FaBook,
  FaCashRegister,
  FaChartBar,
  FaCheckCircle,
  FaClock,
  FaLayerGroup,
  FaPaperPlane,
  FaPencilAlt,
  FaProjectDiagram,
  FaRegCheckCircle,
  FaRegUser,
  FaSearch,
  FaSpinner,
  FaTasks,
  FaUsers,
} from "react-icons/fa";
import { GiMaterialsScience } from "react-icons/gi";
import { MdOutlineLayers, MdOutlineLock } from "react-icons/md";

import { SECTION_TYPES } from "./section-types";
import { STAGES } from "./stages";

export const TITLE_ICON_COMPONENTS = {
  [SECTION_TYPES.LiteratureReview]: FaBook,
  [SECTION_TYPES.Plan]: FaTasks,
  [SECTION_TYPES.Report]: FaChartBar,
  [SECTION_TYPES.ProjectGroup]: FaProjectDiagram,
  [SECTION_TYPES.Note]: GiMaterialsScience,
  Project: FaLayerGroup,
  Students: FaUsers,
  Instructors: FaRegUser,
  ProjectType: MdOutlineLayers,
  RegistrationRule: FaCashRegister,
};

export const STATUS_ICON_COMPONENTS = {
  [STAGES.Draft]: { icon: FaPencilAlt, color: "gray" },
  [STAGES.InReview]: { icon: FaSearch, color: "purple" },
  [STAGES.AwaitingChanges]: { icon: FaClock, color: "orange" },
  [STAGES.Approved]: { icon: FaCheckCircle, color: "green" },
  [STAGES.OnGoing]: { icon: FaSpinner, color: "blue.500" },
  [STAGES.Submitted]: { icon: FaRegCheckCircle, color: "green" },
  [STAGES.Locked]: { icon: MdOutlineLock, color: "yellow.500" },
  [STAGES.FeedbackRequested]: { icon: FaPaperPlane, color: "purple" },
  [STAGES.InProgress]: { icon: FaSpinner, color: "blue.700" },
  [STAGES.InProgressPostFeedback]: { icon: FaSpinner, color: "blue.800" },
  [STAGES.Ready]: { icon: FaCheckCircle, color: "green" },
  [STAGES.Deprecated]: { icon: MdOutlineLock, color: "red" },
};

export const TOAST_DEFAULTS = {
  position: "top",
  duration: 2000,
  isClosable: true,
};
