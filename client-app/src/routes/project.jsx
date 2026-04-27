import { useSectionsListByProject } from "api";
import {
  EXPERIMENTS_PERMISSIONS,
  PROJECTMANAGEMENT_PERMISSIONS,
  SECTION_TYPES as EXPERIMENT_DATA_TYPES,
} from "constants";
import { useUser } from "contexts";
import { ProtectedRoutes } from "layouts/protected-routes";
import { NotFound } from "pages/error";
// import {
//   LiteratureReviewOverview,
//   NoteOverview,
//   PlanOverview,
//   ReportOverview,
// } from "pages/experiment/overview";
// import {
//   GroupProjectSummarySection,
//   LiteratureReviewSection,
//   NoteSection,
//   PlanSection,
//   ReportSection,
// } from "pages/experiment/section";
// import { StudentExperimentList } from "pages/experiment/summary";
import {
  // ProjectGroupStudentList,
  ProjectInstructorList,
  ProjectList,
  ProjectOverview,
} from "pages/project";
import { useEffect } from "react";
import { Route, Routes, useNavigate, useParams } from "react-router-dom";

export const Project = () => {
  const { user } = useUser();

  return (
    <Routes>
      <Route path="/" element={<ProjectList />} />
      <Route
        path=":projectId"
        element={
          <ProtectedRoutes
            isAuthorized={(user) =>
              [
                EXPERIMENTS_PERMISSIONS.ViewExperiments,
                EXPERIMENTS_PERMISSIONS.ViewProjectExperiments,
              ].some((permission) => user.permissions?.includes(permission))
            }
          />
        }
      >
        <Route
          index
          element={
            [
              PROJECTMANAGEMENT_PERMISSIONS.CreateProjects,
              PROJECTMANAGEMENT_PERMISSIONS.CreateProjectGroups,
            ].some((permission) => user.permissions?.includes(permission)) ? (
              <ProjectOverview />
            ) : (
              <StudentExperimentList />
            )
          }
        />
      </Route>

      <Route
        path=":projectId/project-groups/:projectGroupId/students/:studentId"
        element={
          <ProtectedRoutes
            isAuthorized={(user) =>
              [
                EXPERIMENTS_PERMISSIONS.ViewProjectGroupExperiments,
                EXPERIMENTS_PERMISSIONS.ViewProjectExperiments,
              ].some((permission) => user.permissions?.includes(permission))
            }
          />
        }
      >
        <Route index element={<StudentExperimentList />} />
      </Route>

      <Route
        path=":projectId/project-groups/:projectGroupId/literature-reviews/:literatureReviewId/overview"
        element={
          <ProtectedRoutes
            isAuthorized={(user) =>
              [
                EXPERIMENTS_PERMISSIONS.ViewExperiments,
                EXPERIMENTS_PERMISSIONS.ViewProjectExperiments,
              ].some((permission) => user.permissions?.includes(permission))
            }
          />
        }
      >
        <Route
          index
          element={
            <>
              <RedirectToSectionForm
                sectionType={EXPERIMENT_DATA_TYPES.LiteratureReview}
              />
              <LiteratureReviewOverview />
            </>
          }
        />
      </Route>

      <Route
        path=":projectId/project-groups/:projectGroupId/plans/:planId/overview"
        element={
          <ProtectedRoutes
            isAuthorized={(user) =>
              [
                EXPERIMENTS_PERMISSIONS.ViewExperiments,
                EXPERIMENTS_PERMISSIONS.ViewProjectExperiments,
              ].some((permission) => user.permissions?.includes(permission))
            }
          />
        }
      >
        <Route index element={<PlanOverview />} />
      </Route>

      <Route
        path=":projectId/project-groups/:projectGroupId/notes/:noteId/overview"
        element={
          <ProtectedRoutes
            isAuthorized={(user) =>
              [
                EXPERIMENTS_PERMISSIONS.ViewExperiments,
                EXPERIMENTS_PERMISSIONS.ViewProjectExperiments,
              ].some((permission) => user.permissions?.includes(permission))
            }
          />
        }
      >
        <Route index element={<NoteOverview />} />
      </Route>

      <Route
        path=":projectId/project-groups/:projectGroupId/reports/:reportId/overview"
        element={
          <ProtectedRoutes
            isAuthorized={(user) =>
              [
                EXPERIMENTS_PERMISSIONS.ViewExperiments,
                EXPERIMENTS_PERMISSIONS.ViewProjectExperiments,
              ].some((permission) => user.permissions?.includes(permission))
            }
          />
        }
      >
        <Route index element={<ReportOverview />} />
      </Route>

      <Route
        path=":projectId/project-groups/:projectGroupId/literature-reviews/:literatureReviewId/sections/:sectionId"
        element={
          <ProtectedRoutes
            isAuthorized={(user) =>
              [
                EXPERIMENTS_PERMISSIONS.ViewExperiments,
                EXPERIMENTS_PERMISSIONS.ViewProjectExperiments,
              ].some((permission) => user.permissions?.includes(permission))
            }
          />
        }
      >
        <Route index element={<LiteratureReviewSection />} />
      </Route>

      <Route
        path=":projectId/project-groups/:projectGroupId/plans/:planId/sections/:sectionId"
        element={
          <ProtectedRoutes
            isAuthorized={(user) =>
              [
                EXPERIMENTS_PERMISSIONS.ViewExperiments,
                EXPERIMENTS_PERMISSIONS.ViewProjectExperiments,
              ].some((permission) => user.permissions?.includes(permission))
            }
          />
        }
      >
        <Route index element={<PlanSection />} />
      </Route>

      <Route
        path=":projectId/project-groups/:projectGroupId/reports/:reportId/sections/:sectionId"
        element={
          <ProtectedRoutes
            isAuthorized={(user) =>
              [
                EXPERIMENTS_PERMISSIONS.ViewExperiments,
                EXPERIMENTS_PERMISSIONS.ViewProjectExperiments,
              ].some((permission) => user.permissions?.includes(permission))
            }
          />
        }
      >
        <Route index element={<ReportSection />} />
      </Route>

      <Route
        path=":projectId/project-groups/:projectGroupId/notes/:noteId/sections/:sectionId"
        element={
          <ProtectedRoutes
            isAuthorized={(user) =>
              [
                EXPERIMENTS_PERMISSIONS.ViewExperiments,
                EXPERIMENTS_PERMISSIONS.ViewProjectExperiments,
              ].some((permission) => user.permissions?.includes(permission))
            }
          />
        }
      >
        <Route index element={<NoteSection />} />
      </Route>
      <Route
        path=":projectId/project-groups/:projectGroupId/students"
        element={
          <ProtectedRoutes
            isAuthorized={(user) =>
              user.permissions?.includes(
                EXPERIMENTS_PERMISSIONS.ViewProjectGroupExperiments,
              )
            }
          />
        }
      >
        <Route index element={<ProjectGroupStudentList />} />
      </Route>

      <Route
        path=":projectId/project-groups/:projectGroupId/activities"
        element={
          <ProtectedRoutes
            isAuthorized={(user) =>
              [
                EXPERIMENTS_PERMISSIONS.ViewExperiments,
                EXPERIMENTS_PERMISSIONS.ViewProjectExperiments,
              ].some((permission) => user.permissions?.includes(permission))
            }
          />
        }
      >
        <Route index element={<GroupProjectSummarySection />} />
      </Route>

      <Route
        path=":projectId/instructors"
        element={
          <ProtectedRoutes
            isAuthorized={(user) =>
              user.permissions?.includes(
                PROJECTMANAGEMENT_PERMISSIONS.InviteInstructors,
              )
            }
          />
        }
      >
        <Route index element={<ProjectInstructorList />} />
      </Route>

      <Route path="*" element={<NotFound />} />
    </Routes>
  );
};

/**
 * Get the path to the overview page
 */
export const buildOverviewPath = (
  sectionType,
  projectId,
  projectGroupId,
  recordId,
) => {
  const segment = sectionTypePathMap[sectionType];
  return `/projects/${projectId}/project-groups/${projectGroupId}/${segment}/${recordId}/overview`;
};

/**
 * Get the path to the section page
 */
export const buildSectionFormPath = (
  sectionType,
  projectId,
  projectGroupId,
  recordId,
  sectionId,
) => {
  const segment = sectionTypePathMap[sectionType];
  return `/projects/${projectId}/project-groups/${projectGroupId}/${segment}/${recordId}/sections/${sectionId}`;
};

/**
 * Get the path to the project page
 */
export const buildProjectPath = (
  projectId,
  projectGroupId,
  studentId,
  isInstructorOrGroupMember = projectGroupId && studentId,
) => {
  return isInstructorOrGroupMember
    ? `/projects/${projectId}/project-groups/${projectGroupId}/students/${studentId}`
    : `/projects/${projectId}`;
};

/**
 * Get the path to the activities page
 */
export const buildActivitiesPath = (projectId, projectGroupId) => {
  return `/projects/${projectId}/project-groups/${projectGroupId}/activities`;
};

/**
 * Get the path to the students project group page
 */
export const buildStudentsProjectGroupPath = (projectId, projectGroupId) => {
  return `/projects/${projectId}/project-groups/${projectGroupId}/students`;
};

/**
 * Redirect to the section form route if there is only one section for the given section type.
 * Note: for now only applied to the literature review overview route but can be used in other routes if needed.
 */
const RedirectToSectionForm = ({ sectionType }) => {
  const navigate = useNavigate();
  const { projectId, projectGroupId } = useParams();
  const recordId = useParams()[sectionTypeToIdMap[sectionType]];

  const { data: projectSections } = useSectionsListByProject(projectId);

  useEffect(() => {
    if (!projectSections || !sectionType) return;

    const matchingSections = projectSections.filter(
      (section) => section.sectionType.name === sectionType,
    );

    if (matchingSections.length === 1) {
      const nextPath = buildSectionFormPath(
        sectionType,
        projectId,
        projectGroupId,
        recordId,
        matchingSections[0].id,
      );
      navigate(nextPath, { replace: true });
    }
  }, [
    navigate,
    projectGroupId,
    projectId,
    projectSections,
    recordId,
    sectionType,
  ]);

  return null;
};

const sectionTypeToIdMap = {
  LiteratureReview: "literatureReviewId",
  Plan: "planId",
  Note: "noteId",
  Report: "reportId",
};

const sectionTypePathMap = {
  [EXPERIMENT_DATA_TYPES.LiteratureReview]: "literature-reviews",
  [EXPERIMENT_DATA_TYPES.Plan]: "plans",
  [EXPERIMENT_DATA_TYPES.Note]: "notes",
  [EXPERIMENT_DATA_TYPES.Report]: "reports",
};
