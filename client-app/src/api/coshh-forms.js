import { useBackendApi } from "contexts";
import useSWR from "swr";

export const fetchKeys = {
  coshhFormsList: "templates",
  coshhForm: (coshhFormId) => `templates/${coshhFormId}`,
  sectionTypesList: "templates/section_types",
};

export const getCoshhFormsApi = ({ api }) => ({
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

export const useCoshhFormsList = () => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    fetchKeys.coshhFormsList,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    {
      suspense: true,
    },
  );
};

export const useCoshhForm = (coshhFormId) => {
  const { apiFetcher } = useBackendApi();
  return useSWR(
    coshhFormId ? fetchKeys.coshhForm(coshhFormId) : null,
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
