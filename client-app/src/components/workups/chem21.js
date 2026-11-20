/*
 * CHEM21 solvent flags, as stored in the Solvent table, and how they're shown.
 * The colours are each user's own (they can change them in their account), so
 * they come from the backend with the solvent list. These match the reaction
 * table's colours (static/js/custom_colours.js).
 */

const FLAGS = {
  1: { label: "Highly hazardous", colourKey: "HighlyHazardous" },
  2: { label: "Hazardous", colourKey: "Hazardous" },
  3: { label: "Problematic", colourKey: "Problematic" },
  4: { label: "Recommended", colourKey: "Recommended" },
  5: { label: "Non-CHEM21", colourKey: null },
};

const NON_CHEM21_STYLE = { bg: "white", color: "black" };

/**
 * The label for a CHEM21 flag.
 * @param {number | null} flag The flag, 1 (highly hazardous) to 5 (not in CHEM21).
 * @returns {string} The label, or "Non-CHEM21" for an unknown flag.
 */
export const flagLabel = (flag) => FLAGS[flag]?.label ?? FLAGS[5].label;

/**
 * The background and text colours for a CHEM21 flag.
 * @param {number | null} flag The flag.
 * @param {object | undefined} colours The user's colours, e.g.
 *   `{ Recommended: "#00ff00", Recommended_text: "#000000", ... }`.
 * @returns {{ bg: string, color: string }} Chakra style props.
 */
export const flagStyle = (flag, colours) => {
  const colourKey = FLAGS[flag]?.colourKey;
  if (!colourKey || !colours?.[colourKey]) return NON_CHEM21_STYLE;
  return {
    bg: colours[colourKey],
    color: colours[`${colourKey}_text`] ?? "black",
  };
};
