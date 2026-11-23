/*
 * Mixtures, e.g. a 90:10 hexane/EtOAc eluent. Each material in the mixture has
 * a percentage, and the user enters one total volume for the whole mixture,
 * which is split into each material's own amount (so green chemistry metrics
 * still get a volume per solvent). Values are kept as typed strings, like the
 * rest of the step form.
 */

// what a mixture's percentages add up to (matches the backend)
export const PROPORTION_TOTAL = 100;

/** Rounds a percentage to 2 decimal places, as the backend stores it. */
const roundPercentage = (value) => Math.round(value * 100) / 100;

/** Rounds an amount to 4 decimal places, as the backend stores it. */
const roundAmount = (value) => Math.round(value * 10000) / 10000;

/**
 * Whether a material row has a percentage entered.
 * @param {object} row A material row.
 * @returns {boolean}
 */
export const hasProportion = (row) =>
  row.proportion !== "" &&
  row.proportion !== null &&
  row.proportion !== undefined;

/**
 * The sum of the rows' entered percentages.
 * @param {object[]} rows Material rows.
 * @returns {number}
 */
export const proportionTotal = (rows) =>
  roundPercentage(
    rows.reduce((sum, row) => sum + (Number(row.proportion) || 0), 0),
  );

/**
 * Whether every row has a percentage and they add up to 100, so a volume can
 * be split between them.
 * @param {object[]} rows Material rows.
 * @returns {boolean}
 */
export const isProportionComplete = (rows) =>
  rows.every(hasProportion) && proportionTotal(rows) === PROPORTION_TOTAL;

/**
 * Keeps a mixture's percentages adding up to 100 by setting one row to whatever
 * is left: the last row, or the one before it when the last is being edited
 * (so with two materials, typing 90 in either sets the other to 10). A single
 * material has no percentage.
 * @param {object[]} rows The slot's material rows, in order.
 * @param {string | null} editedKey The row the user just changed, if any.
 * @returns {object[]} The rows, with the balancing row updated.
 */
export const balanceProportions = (rows, editedKey) => {
  if (rows.length < 2) return rows.map((row) => ({ ...row, proportion: "" }));
  const lastIndex = rows.length - 1;
  const balanceIndex =
    rows[lastIndex].key === editedKey ? lastIndex - 1 : lastIndex;
  const others = proportionTotal(
    rows.filter((_, index) => index !== balanceIndex),
  );
  const remainder = roundPercentage(PROPORTION_TOTAL - others);
  // leave it to the user if the others already add up to more than 100
  if (remainder < 0) return rows;
  return rows.map((row, index) =>
    index === balanceIndex ? { ...row, proportion: String(remainder) } : row,
  );
};

/**
 * The mixture's total volume, worked out from its materials' amounts (e.g.
 * when a saved step is opened), or "" if none are known.
 * @param {object[]} rows The mixture's material rows.
 * @returns {string}
 */
export const mixtureVolume = (rows) => {
  const amounts = rows
    .map((row) => row.quantity)
    .filter((quantity) => quantity !== "" && quantity !== null);
  if (amounts.length === 0) return "";
  return String(
    roundAmount(amounts.reduce((sum, amount) => sum + Number(amount), 0)),
  );
};

/**
 * Splits a mixture's total volume between its materials by their
 * percentages. Leaves the amounts alone while the percentages are incomplete,
 * and clears them when the volume is cleared.
 * @param {object[]} rows The mixture's material rows.
 * @param {string} volume The total volume, as typed.
 * @returns {object[]} The rows with their amounts set.
 */
export const splitVolume = (rows, volume) => {
  if (volume === "") return rows.map((row) => ({ ...row, quantity: "" }));
  const total = Number(volume);
  if (!Number.isFinite(total) || total < 0 || !isProportionComplete(rows)) {
    return rows;
  }
  return rows.map((row) => ({
    ...row,
    quantity: String(
      roundAmount((total * Number(row.proportion)) / PROPORTION_TOTAL),
    ),
  }));
};
