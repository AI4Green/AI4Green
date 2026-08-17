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
    {
      suspense: true,
      fallbackData: [
        { id: "1", name: "Example Type A" },
        { id: "2", name: "Example Type B" },
      ],
    },
  );
};

export const useProjectType = (projectTypeId) => {
  const { apiFetcher } = useBackendApi();
  const mockData = {
    id: "mock-id-123",
    name: "Green Chemistry Research",
    description:
      "A specialized project type for analyzing sustainable chemical reactions and waste reduction metrics.",
    // Add any other fields your 'Section' or 'Area' components expect
    sections: [],
  };

  return useSWR(
    projectTypeId ? fetchKeys.projectType(projectTypeId) : null,
    async (url) => {
      // const data = await apiFetcher(url);

      return mockData;
    },
    { suspense: true },
  );
};

export const useSectionTypesList = () => {
  const { apiFetcher } = useBackendApi();

  const mockSectionTypes = [
    {
      id: "section-type-1",
      name: "Standard Analysis",
      description: "Basic green chemistry metrics",
      order: 1,
    },
    {
      id: "section-type-2",
      name: "Solvent Selection",
      description: "Safety and environmental impact of solvents",
      order: 2,
    },
  ];

  return useSWR(
    fetchKeys.sectionTypesList,
    async (url) => {
      // const data = await apiFetcher(url);
      return mockSectionTypes;
    },
    { suspense: true },
  );
};
