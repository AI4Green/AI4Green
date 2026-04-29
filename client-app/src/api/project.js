import { useBackendApi } from "contexts";
import useSWR from "swr";

export const fetchKeys = {
  list: "template-instances/",
  get: (id) => `template-instances/${id}`,
  studentSummary: (id, studentId) =>
    `projects/${id}/summary${studentId ? `?studentId=${studentId}` : ""}`,
  instructors: (id) => `projects/${id}/instructors`,
  validateInstructor: (id) => `projects/${id}/validate-instructor`,
};

export const getProjectsApi = ({ api }) => ({
  create: (values) =>
    api.post("template-instances/", {
      json: values,
    }),

  edit: (id, values) =>
    api.put(`template-instances/${id}`, {
      json: values,
    }),

  delete: (id) => api.delete(`template-instances/${id}`),

  inviteInstructors: (id, { emails }) =>
    api.post(`template-instances/${id}/invite-instructors`, {
      json: { emails },
    }),

  removeInstructor: (id, instructorId) =>
    api.post(`template-instances/${id}/remove-instructor`, {
      json: { id: instructorId },
    }),
});

export const useProjectsList = () => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    fetchKeys.list,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};

export const useProject = (id) => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    id ? fetchKeys.get(id) : null,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};

export const useProjectSummaryByStudent = (id, studentId) => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    id ? fetchKeys.studentSummary(id, studentId) : null,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};

export const useProjectInstructors = (id) => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    id ? fetchKeys.instructors(id) : null,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};

export const useIsProjectInstructor = (id) => {
  const { api } = useBackendApi();

  return useSWR(
    id ? fetchKeys.validateInstructor(id) : null,
    async () => {
      try {
        await api.post(fetchKeys.validateInstructor(id));
        return true;
      } catch {
        return false;
      }
    },
    { suspense: true },
  );
};
