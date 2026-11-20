/**
 * Picks the definitions for the given keys, in the given order, skipping any
 * the backend doesn't define.
 * @param {object[]} fields The step type's field definitions.
 * @param {string[]} keys The keys wanted.
 * @returns {object[]} The matching definitions.
 */
export const pickFields = (fields, keys) =>
  keys.map((key) => fields.find((field) => field.key === key)).filter(Boolean);
