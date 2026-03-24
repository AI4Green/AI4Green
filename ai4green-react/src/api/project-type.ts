import { useBackendApi } from "contexts";
import useSWR from "swr";

export const fetchKeys = {
  projectTypesList: "project-types/",
  projectType: (projectTypeId) => `project-types/${projectTypeId}`,
  sectionTypesList: "project-types/section-types",
};

export const getProjectTypesApi = ({ api }) => ({
  create: ({ values }) =>
    api.post("project-types/", {
      json: values,
    }),

  edit: ({ values, id }) =>
    api.put(`project-types/${id}`, {
      json: values,
    }),

  delete: (id) => api.delete(`project-types/${id}`),

  advanceStage: (id, stageName) =>
    api.post(`project-types/${id}/advance`, { json: { stageName } }),
});

export const useProjectTypesList = () => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    fetchKeys.projectTypesList,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};

export const useProjectType = (projectTypeId) => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    projectTypeId ? fetchKeys.projectType(projectTypeId) : null,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};

export const useSectionTypesList = () => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    fetchKeys.sectionTypesList,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};
