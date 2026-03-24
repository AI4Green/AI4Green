import { useBackendApi } from "contexts";
import useSWR from "swr";

export const fetchKeys = {
  fieldByName: (projectId, sectionType, fieldName) =>
    `fields/${projectId}/${sectionType}/${fieldName}`,

  sectionFields: (id) => `fields/section/${id}`,

  inputTypes: "input-types",
};

export const getFieldsApi = ({ apiFetcher, api }) => ({
  getFieldByName: async (projectId, sectionType, fieldName) =>
    apiFetcher(fetchKeys.fieldByName(projectId, sectionType, fieldName)),

  save: (id, values) => api.post(`fields/${id}/save`, { json: values }),
});

export const useSectionFields = (id) => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    id ? fetchKeys.sectionFields(id) : null,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};

export const useInputTypes = () => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    fetchKeys.inputTypes,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};
