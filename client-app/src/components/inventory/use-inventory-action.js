import { useToast } from "@chakra-ui/react";
import { useWorkgroupInventory } from "api";
import { GLOBAL_PARAMETERS } from "constants";
import { useState } from "react";

// longest server validation message shown to the user
const MAX_DETAIL_LENGTH = 300;

/**
 * Reads the `message` from a JSON error response. Only a short plain-string
 * message is used, so an HTML error page or stack trace is never shown.
 * @param {Response} response The failed response.
 * @returns {Promise<string | null>} The message, or null if there isn't a usable one.
 */
const readErrorMessage = async (response) => {
  try {
    const { message } = await response.json();
    return typeof message === "string" && message.length <= MAX_DETAIL_LENGTH
      ? message
      : null;
  } catch {
    return null;
  }
};

/**
 * Turns a failed API call into a message for the user.
 * @param {Error} error The error thrown by ky.
 * @returns {Promise<string>} The message to show.
 */
export const describeError = async (error) => {
  const status = error?.response?.status;
  if (status === 404) {
    // the backend also returns 404 for rows in other workgroups
    return "This couldn't be found. It may have been removed, or you may no longer have access to it. Try refreshing the page.";
  }
  if (status === 403) {
    return "Only principal investigators and senior researchers can change this workgroup's inventory.";
  }
  if (status === 400) {
    const detail = await readErrorMessage(error.response);
    return detail ?? "Check the details and try again.";
  }
  return "Something went wrong. Please try again.";
};

/**
 * Runs an inventory write for a modal: tracks the busy state, refreshes the
 * workgroup's inventory and shows a success toast when it works, and keeps
 * an error message for the modal to show when it doesn't.
 * @param {string} workgroupName The workgroup the inventory belongs to.
 * @returns {{
 *   run: (request: () => Promise<unknown>, successTitle: string) => Promise<boolean>,
 *   isLoading: boolean,
 *   error: string | null,
 *   clearError: () => void,
 * }} `run` resolves to true if the request succeeded.
 */
export const useInventoryAction = (workgroupName) => {
  const { mutate } = useWorkgroupInventory(workgroupName);
  const toast = useToast();
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);

  const run = async (request, successTitle) => {
    setIsLoading(true);
    setError(null);
    try {
      await request();
      await mutate();
      toast({
        title: successTitle,
        status: "success",
        duration: GLOBAL_PARAMETERS.ToastDuration,
        isClosable: true,
        position: "top",
      });
      return true;
    } catch (e) {
      console.error(e);
      setError(await describeError(e));
      return false;
    } finally {
      setIsLoading(false);
    }
  };

  return { run, isLoading, error, clearError: () => setError(null) };
};
