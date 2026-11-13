import { useBackendApi } from "contexts";
import useSWR from "swr";

/*
 * Search API (Webapp/sources/blueprints/api/search).
 *
 * Only the reactions in workbooks the user belongs to are searched. The
 * search is a POST with a JSON body because a drawn structure can be too long
 * for a URL; it doesn't change anything.
 */

// not exported, so it doesn't clash with the other modules' fetchKeys in "api"
const fetchKeys = {
  workgroups: "search/workgroups",
};

/** The structure search modes the backend accepts. */
export const STRUCTURE_SEARCH_MODES = ["exact", "substructure"];

/** The parts of a reaction a structure can be matched in. */
export const STRUCTURE_SEARCH_ROLES = [
  "reactant",
  "reagent",
  "solvent",
  "product",
];

/**
 * The user's workgroups and the workbooks they belong to, for the search
 * filters: `{ workgroups: [{ name, workbooks: [name] }] }`.
 */
export const useSearchWorkgroups = () => {
  const { apiFetcher } = useBackendApi();
  return useSWR(fetchKeys.workgroups, apiFetcher, { suspense: true });
};

/** Search actions, added to the backend API context as `search`. */
export const getSearchApi = ({ api }) => ({
  /**
   * Search the user's reactions. Every field is optional, but at least one of
   * keyword, reagent, dateFrom, dateTo or structure is needed.
   * @param {{
   *   workgroup?: string, workbook?: string, keyword?: string,
   *   reagent?: string, dateFrom?: string, dateTo?: string,
   *   structure?: { smiles: string, mode: "exact" | "substructure",
   *                 roles: string[] }
   * }} criteria What to search for. Dates are YYYY-MM-DD.
   * @returns {Promise<{ results: object[], truncated: boolean }>} Each result
   *   has `id, reactionId, name, description, workgroup, workbook,
   *   creatorUsername, timeOfCreation, timeOfUpdate, completionStatus, image,
   *   matchedRoles`.
   */
  searchReactions: (criteria) =>
    api.post("search/reactions", { json: criteria }).json(),
});
