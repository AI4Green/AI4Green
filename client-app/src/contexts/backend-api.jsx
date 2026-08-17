import {
  getAccountApi,
  // getCommentsApi,
  getFieldsApi,
  // getLiteratureReviewsApi,
  // getNotesApi,
  // getPlansApi,
  // getPredictionsApi,
  // getProjectGroupsApi,
  getProjectsApi,
  getProjectTypesApi,
  // getReactionTableApi,
  // getRegistrationRulesApi,
  // getReportsApi,
  getSectionsApi,
  getUserApi,
} from "api";
import ky from "ky";
import { createContext, useCallback, useContext, useMemo } from "react";
import { useTranslation } from "react-i18next";

const BackendApiContext = createContext({});

export const useBackendApi = () => useContext(BackendApiContext);

/** Default KY instance options for hitting the backend API */
export const getBackendDefaults = (language) => ({
  prefixUrl: "/api/",
  credentials: "include",
  headers: {
    "Accept-Language": language,
  },
  timeout: 5 * 60 * 1000,
});

export const BackendApiProvider = ({ children }) => {
  // guarantees i18n is initialised
  // and gives us the instance so we can
  // tell the backend the language in use :)
  const { i18n } = useTranslation();

  // preconfigured ky instance for hitting the backend api
  const api = useMemo(
    () => ky.create(getBackendDefaults(i18n.language)),
    [i18n],
  );

  /**
   * A default fetcher for SWR to get data from the backend API
   * @param {*} path the url path relative to `https://{backend}/api/`
   * @returns
   */
  const apiFetcher = useCallback(
    async (path) => await api.get(path).json(),
    [api],
  );

  const baseContext = useMemo(() => ({ api, apiFetcher }), [api, apiFetcher]);

  const context = useMemo(
    () => ({
      ...baseContext,
      account: getAccountApi(baseContext),
      users: getUserApi(baseContext),
      // registrationRules: getRegistrationRulesApi(baseContext),
      projects: getProjectsApi(baseContext),
      // projectGroups: getProjectGroupsApi(baseContext),
      projectTypes: getProjectTypesApi(baseContext),
      // plans: getPlansApi(baseContext),
      // notes: getNotesApi(baseContext),
      // literatureReviews: getLiteratureReviewsApi(baseContext),
      // reports: getReportsApi(baseContext),
      // reactionTable: getReactionTableApi(baseContext),
      // predictions: getPredictionsApi(baseContext),
      // comments: getCommentsApi(baseContext),
      sections: getSectionsApi(baseContext),
      fields: getFieldsApi(baseContext),
    }),
    [baseContext],
  );

  return (
    <BackendApiContext.Provider value={context}>
      {children}
    </BackendApiContext.Provider>
  );
};
