import { SECTION_TYPES } from "./section-types";

export const STAGES = {
  Draft: "Draft",
  InReview: "In Review",
  AwaitingChanges: "Awaiting Changes",
  Submitted: "Submitted",
  Approved: "Approved",
  OnGoing: "On Going",
  Locked: "Locked",
  InProgress: "In Progress",
  FeedbackRequested: "Feedback Requested",
  InProgressPostFeedback: "In Progress Post-Feedback",
  Ready: "Ready",
  Deprecated: "Deprecated",
};

export const STAGE_TYPES = {
  ProjectType: "ProjectType",
  LiteratureReview: SECTION_TYPES.LiteratureReview,
  Plan: SECTION_TYPES.Plan,
  Note: SECTION_TYPES.Note,
  Report: SECTION_TYPES.Report,
};
