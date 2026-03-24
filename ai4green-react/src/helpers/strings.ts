export const capitalise = (s) =>
  s
    .split("")
    .map((x, i) => (i ? x : x.toUpperCase()))
    .join("");

/**
 * Convert a "Yes/No" string to a boolean,
 * where "Yes" is true and everything else is false.
 * @param {*} value
 * @returns
 */
export const radioBool = (value) => value === "Yes";

/**
 * Count the words in a string
 * @param {*} str
 * @returns
 */
export const countWords = (str) => {
  if (str.length === 0) return 0;
  return str.trim().split(/\s+/).length;
};
