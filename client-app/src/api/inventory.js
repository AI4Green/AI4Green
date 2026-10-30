import { useBackendApi } from "contexts";
import useSWR from "swr";

export const fetchKeys = {
  workgroup: (workgroupName) =>
    `inventory/${encodeURIComponent(workgroupName)}`,
};

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
