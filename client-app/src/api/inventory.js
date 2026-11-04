import { useBackendApi } from "contexts";
import useSWR from "swr";

/*
 * Inventory API (Webapp/sources/blueprints/api/inventory).
 *
 * Any workgroup member can read; only principal investigators and senior
 * researchers can write (the backend returns 403 otherwise, and the inventory
 * payload's `permissions.canEdit` tells the UI which actions to show).
 *
 * Every write sends a JSON body. A cross-site HTML form can't, and app-wide
 * CSRF tokens are currently disabled, so the backend rejects anything that
 * isn't `application/json`.
 */

/** Builds the base path for a workgroup's inventory. */
const workgroupPath = (workgroupName) =>
  `inventory/${encodeURIComponent(workgroupName)}`;

/**
 * Checks an id is a positive integer before it goes into a URL, so a bad value
 * (e.g. "../users") can't point the request at a different API route.
 * @param {*} id An item or location id.
 * @returns {number} The id.
 * @throws {TypeError} If the id isn't a positive integer.
 */
const safeId = (id) => {
  if (!Number.isInteger(id) || id <= 0) {
    throw new TypeError(`Invalid inventory id: ${id}`);
  }
  return id;
};

/** Builds the path for one of a workgroup's containers. */
const itemPath = (workgroupName, itemId) =>
  `${workgroupPath(workgroupName)}/items/${safeId(itemId)}`;

/** Builds the path for one of a workgroup's locations. */
const locationPath = (workgroupName, locationId) =>
  `${workgroupPath(workgroupName)}/locations/${safeId(locationId)}`;

export const fetchKeys = {
  workgroup: (workgroupName) => workgroupPath(workgroupName),
  itemHistory: (workgroupName, itemId) =>
    `${itemPath(workgroupName, itemId)}/history`,
};

/**
 * Inventory write actions, added to the backend API context as `inventory`.
 * Every call is scoped to a workgroup, which the backend checks membership of.
 * Quantities are always in the container's own unit.
 */
export const getInventoryApi = ({ api }) => ({
  /**
   * Look up a compound by name or CAS number, in the compound database and the
   * novel compounds of the workgroup's workbooks. Read-only. Resolves to one of:
   * - `{ status: "found", compound }`
   * - `{ status: "cas_not_found", cas }`: a CAS number with no match
   * - `{ status: "not_found", suggestions }`: a name with no exact match, and
   *   up to a few similar compound names
   * `compound` is `{ source: "compound" | "novel", id, name, cas, hazards,
   * workbook, workgroupCompoundId }`; `workbook` is set for novel compounds and
   * `workgroupCompoundId` when the workgroup already holds the compound.
   */
  lookupCompound: (workgroupName, query) =>
    api
      .get(`${workgroupPath(workgroupName)}/compounds/lookup`, {
        searchParams: { query },
      })
      .json(),

  /**
   * Add a novel compound to one of the workgroup's workbooks, as the reaction
   * table does. Resolves to the new compound, in the same shape as a lookup.
   */
  addNovelCompound: (workgroupName, values) =>
    api
      .post(`${workgroupPath(workgroupName)}/novel-compounds`, { json: values })
      .json(),

  /**
   * Add a container. Pass `workgroupCompoundId` for a compound the workgroup
   * already holds, or `compound: { source, id }` from a lookup to start
   * holding a new one.
   */
  addItem: (workgroupName, values) =>
    api.post(`${workgroupPath(workgroupName)}/items`, { json: values }),

  /** Set a container's measured quantity; the backend records the difference. */
  adjustItem: (workgroupName, itemId, { quantity, notes }) =>
    api.post(`${itemPath(workgroupName, itemId)}/adjust`, {
      json: { quantity, notes },
    }),

  moveItem: (workgroupName, itemId, { locationId }) =>
    api.post(`${itemPath(workgroupName, itemId)}/move`, {
      json: { locationId },
    }),

  /** Dispose of a container: empties it and archives it. */
  disposeItem: (workgroupName, itemId, { notes }) =>
    api.post(`${itemPath(workgroupName, itemId)}/dispose`, {
      json: { notes },
    }),

  addLocation: (workgroupName, { name, description }) =>
    api.post(`${workgroupPath(workgroupName)}/locations`, {
      json: { name, description },
    }),

  editLocation: (workgroupName, locationId, { name, description }) =>
    api.put(locationPath(workgroupName, locationId), {
      json: { name, description },
    }),

  setLocationArchived: (workgroupName, locationId, isArchived) =>
    api.patch(locationPath(workgroupName, locationId), {
      json: { isArchived },
    }),
});

/**
 * Get a workgroup's inventory: its storage locations and the compounds it
 * holds, each with its containers under `items`.
 * @param {string} workgroupName The workgroup name
 * @returns
 */
export const useWorkgroupInventory = (workgroupName) => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    workgroupName ? fetchKeys.workgroup(workgroupName) : null,
    async (url) => {
      const data = await apiFetcher(url);
      return data;
    },
    { suspense: true },
  );
};

/**
 * Get a container's ledger: every quantity change, oldest first. Doesn't use
 * suspense, so the history drawer can show its own loading and error states.
 * @param {string} workgroupName The workgroup name
 * @param {number | null} itemId The container id, or null to skip fetching
 * @returns
 */
export const useItemHistory = (workgroupName, itemId) => {
  const { apiFetcher } = useBackendApi();

  return useSWR(
    workgroupName && itemId
      ? fetchKeys.itemHistory(workgroupName, itemId)
      : null,
    apiFetcher,
    // the drawer shows the error; retrying in the background won't fix a 4xx
    { shouldRetryOnError: false },
  );
};
