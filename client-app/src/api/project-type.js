import { useBackendApi } from "contexts";
import useSWR from "swr";

export const fetchKeys = {
  projectTypesList: "coshh/templates",
  projectType: (projectTypeId) => `coshh/templates/${projectTypeId}`,
  sectionTypesList: "coshh/templates/section_types",
};

export const getProjectTypesApi = ({ api }) => ({
  create: ({ values }) =>
    api.post("coshh/templates", {
      json: values,
    }),

  edit: ({ values, id }) =>
    api.put(`coshh/templates/${id}`, {
      json: values,
    }),

  delete: (id) => api.delete(`coshh/templates/${id}`),

  advanceStage: (id, stageName) =>
    api.post(`coshh/templates/${id}/advance`, { json: { stageName } }),
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
