import { array, object, string } from "yup";

export const validationSchema = (projects, projectTypes) =>
  object().shape({
    name: string()
      .notOneOf(
        // fails if project name already exists
        projects.map((project) => project.name),
        "Project name already exist",
      )
      .required("Project name required"),
    projectTypeId: array()
      .length(1, "Please select one project type")
      .of(
        string().oneOf(
          projectTypes.map((projectType) => String(projectType.id)),
          "Invalid project type",
        ),
      )
      .required("Project type is required"),
  });
