import { useBackendApi } from "contexts";
import useSWR from "swr";

const fetchKeys = {
  me: "users/me",
  userList: "users",
  user: (id) => `users/${id}`,
};

/**
 * Get an authenticated user's profile
 * @returns
 */
export const useProfile = () => {
  // const { apiFetcher } = useBackendApi();

  return {
    data: {
      fullName: "Dev Admin",
      email: "admin@dev.local",
      uiCulture: "en",
      // IMPORTANT: Add all permissions you want to test in the UI
      permissions: [
        "CreateProjectTypes",
        "EditProjectTypes",
        "DeleteProjectTypes",
        "ViewProjectTypes",
      ],
    },
    mutate: () => Promise.resolve(),
    isLoading: false,
    isValidating: false,
    error: null,
  };

  // return useSWR(
  //   fetchKeys.me,
  //   async (url) => {
  //     try {
  //       return await apiFetcher(url);
  //     } catch (e) {
  //       if (e?.response?.status === 401) return null;
  //       throw e;
  //     }
  //   },
  //   { suspense: true, refreshInterval: 600000 }, // re-check profile every 10mins in case of cookie expiry
  // );
};

export const getUserApi = ({ api }) => ({
  /**
   * Save the User's UI Culture
   * @param {*} culture a culture/locale/language code
   * @returns
   */
  setUICulture: (culture) =>
    api.put("users/ui-culture", {
      json: culture,
    }),

  // Update user roles
  setUserRoles: ({ id, values }) =>
    api.put(`users/${id}/roles/update`, {
      json: values.roles,
    }),

  // Update user email
  setUserEmail: ({ id, values }) =>
    api.put(`users/${id}/email/request-change?newEmail=${values.email}`),

  // Delete user and email delete update if true
  delete: (user) =>
    api.delete(`users/${user.id}?notify=${user.sendUpdateEmail}`),

  /**
   * Invite a new User
   * @param {*} body
   */
  invite: (values) =>
    api.post("users/invite", {
      json: values,
    }),

  /**
   * Resend an account confirmation email
   * @param {*} userIdOrEmail The User ID or Email Address
   * @returns
   */
  resendInvite: (userIdOrEmail) =>
    api.put("users/invite/resend", {
      json: userIdOrEmail,
    }),
});

/**
 * Get a list of users
 * @returns
 */
export const useUserList = () => {
  const { apiFetcher } = useBackendApi();
  return useSWR(
    fetchKeys.userList,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};

/**
 * Get a user based on their user id
 * @returns
 */
export const useUserSelected = (id) => {
  const { apiFetcher } = useBackendApi();
  return useSWR(
    fetchKeys.user(id),
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};
