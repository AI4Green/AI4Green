import { useBackendApi } from "contexts";
import useSWR from "swr";

export const fetchKeys = {
  fieldByName: (projectId, sectionType, fieldName) =>
    `fields/${projectId}/${sectionType}/${fieldName}`,

  sectionFields: (id) => `sections/${id}/fields`,

  inputTypes: "input-types",
};

export const getFieldsApi = ({ apiFetcher, api }) => ({
  getFieldByName: async (projectId, sectionType, fieldName) =>
    apiFetcher(fetchKeys.fieldByName(projectId, sectionType, fieldName)),

  save: (id, values) => api.post(`sections/${id}/fields`, { json: values }),
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
