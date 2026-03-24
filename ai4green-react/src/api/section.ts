import { useBackendApi } from "contexts";
import useSWR from "swr";

export const fetchKeys = {
  listByProjectType: (id) => `sections/project-type/${id}`,

  listByProject: (id) => `sections/project/${id}`,

  listBySectionType: (id, name) =>
    `sections/section-type/${name}/project/${id}`,

  file: (sectionId, recordId, fileLocation, fileName) =>
    `sections/file?sectionId=${sectionId}&recordId=${recordId}&fileLocation=${fileLocation}&name=${fileName}`,
};

export const getSectionsApi = ({ api }) => ({
  save: (values) => api.post(`sections/save`, { json: values }),

  saveFieldResponses: (formValues) =>
    api.put(`sections/SaveSection`, { body: formValues }),

  downloadSectionFile: (sectionId, recordId, fileLocation, fileName) =>
    api.get(fetchKeys.file(sectionId, recordId, fileLocation, fileName)),
});

export const useSectionsListByProjectType = (id) => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    id ? fetchKeys.listByProjectType(id) : null,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};

export const useSectionsListByProject = (id) => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    id ? fetchKeys.listByProject(id) : null,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};

export const useSectionsListBySectionType = (id, name) => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    id && name ? fetchKeys.listBySectionType(id, name) : null,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};
