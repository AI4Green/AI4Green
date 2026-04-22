import { useBackendApi } from "contexts";
import useSWR from "swr";

export const fetchKeys = {
  projectTypesList: "templates",
  projectType: (projectTypeId) => `templates/${projectTypeId}`,
  sectionTypesList: "templates/section_types",
};

export const getProjectTypesApi = ({ api }) => ({
  create: ({ values }) =>
    api.post("templates", {
      json: values,
    }),

  edit: ({ values, id }) =>
    api.put(`templates/${id}`, {
      json: values,
    }),

  delete: (id) => api.delete(`templates/${id}`),

  advanceStage: (id, stageName) =>
    api.post(`templates/${id}/advance`, { json: { stageName } }),
});

export const useProjectTypesList = () => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    fetchKeys.projectTypesList,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    {
      suspense: true,
    },
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
