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
 * Turns a failed workup API call into a message for the user.
 * @param {Error} error The error thrown by ky.
 * @returns {Promise<string>} The message to show.
 */
export const describeWorkupError = async (error) => {
  const status = error?.response?.status;
  if (status === 400 || status === 403) {
    const detail = await readErrorMessage(error.response);
    if (detail) return detail;
  }
  if (status === 403) {
    return "You can't change this workup.";
  }
  if (status === 404) {
    return "This couldn't be found. It may have been removed. Try refreshing the page.";
  }
  return "Something went wrong. Please try again.";
};
