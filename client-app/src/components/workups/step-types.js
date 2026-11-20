import {
  FaBalanceScale,
  FaClock,
  FaCompressArrowsAlt,
  FaCube,
  FaEyeDropper,
  FaFilter,
  FaFire,
  FaFlask,
  FaGem,
  FaLayerGroup,
  FaMagnet,
  FaPencilAlt,
  FaSync,
  FaThermometerHalf,
  FaTint,
  FaVial,
  FaVials,
  FaWind,
} from "react-icons/fa";

import { ChromatographyStep, GenericStep } from "./step-layouts";

/*
 * How each workup step type looks, keyed by step type (the backend's
 * models.WorkupStepType). This mirrors the COSHH canvas, where each input type
 * has its own appearance. Each entry has:
 * - label, icon, color: how the step is shown in lists and card headers.
 * - Component (optional): lays out the step's inputs. Types without one use
 *   GenericStep, which shows every material slot and setting the backend
 *   defines, so new ones appear without front-end changes.
 * - summary({ values, materials, temperature }) (optional): a one-line
 *   description for a collapsed card. Types without one list their materials
 *   and temperature.
 */

/** Joins the non-empty parts of a summary. */
const joinParts = (parts, separator = " · ") =>
  parts
    .filter((part) => part !== undefined && part !== null && part !== "")
    .join(separator);

/** Formats a value with its unit, or nothing if the value is empty. */
const withUnit = (value, unit) =>
  value === undefined || value === null || value === ""
    ? ""
    : `${value}${unit ? ` ${unit}` : ""}`;

/** Describes the named materials in one slot (or all slots), e.g. "Brine 20 mL". */
const materialsText = (materials, slot) =>
  joinParts(
    materials
      .filter((row) => row.name && (!slot || row.slot === slot))
      .map((row) =>
        joinParts([row.name, withUnit(row.quantity, row.unit)], " "),
      ),
    " + ",
  );

/** "3 ×" in front of a repeated slot's materials, when there are repeats. */
const repeated = (repeats, text) =>
  text && Number(repeats) > 1 ? `${repeats} × ${text}` : text;

/** The default summary: what was used and at what temperature. */
const defaultSummary = ({ materials, temperature }) =>
  joinParts([materialsText(materials), withUnit(temperature, "°C")]);

export const STEP_TYPES_MAP = {
  CUSTOM: {
    label: "Custom step",
    icon: FaPencilAlt,
    color: "gray",
    summary: ({ values, materials }) =>
      joinParts([values.title, materialsText(materials)]),
  },
  ADDITION: {
    label: "Addition",
    icon: FaFlask,
    color: "blue",
    summary: ({ values, materials, temperature }) =>
      joinParts([
        materialsText(materials, "material"),
        values.rate,
        withUnit(temperature, "°C"),
      ]),
  },
  ALIQUOT: {
    label: "Aliquot",
    icon: FaEyeDropper,
    color: "blue",
    summary: ({ values }) =>
      joinParts([
        joinParts([values.amount, values.amount_unit], " "),
        values.purpose,
      ]),
  },
  TEMPERATURE: {
    label: "Temperature change",
    icon: FaThermometerHalf,
    color: "orange",
    summary: ({ values, temperature }) =>
      joinParts([withUnit(temperature, "°C"), values.duration]),
  },
  CONCENTRATION: {
    label: "Concentration",
    icon: FaCompressArrowsAlt,
    color: "purple",
    summary: ({ values, temperature }) =>
      joinParts([
        values.method,
        withUnit(temperature, "°C"),
        withUnit(values.pressure, "mbar"),
      ]),
  },
  EXTRACTION: {
    label: "Extraction",
    icon: FaVials,
    color: "teal",
    summary: ({ values, materials }) =>
      joinParts([
        repeated(values.repeats, materialsText(materials, "solvent")),
        values.phase_kept && `keep ${values.phase_kept.toLowerCase()}`,
      ]),
  },
  FILTRATION: {
    label: "Filtration",
    icon: FaFilter,
    color: "teal",
    summary: ({ values, materials }) =>
      joinParts([
        values.filter_medium,
        materialsText(materials),
        values.phase_kept && `keep ${values.phase_kept.toLowerCase()}`,
      ]),
  },
  WASH: {
    label: "Wash",
    icon: FaTint,
    color: "cyan",
    summary: ({ values, materials }) =>
      joinParts(
        [
          materialsText(materials, "organic_phase"),
          repeated(values.repeats, materialsText(materials, "aqueous_phase")),
        ],
        " / ",
      ),
  },
  DRY_IN_VACUUM: {
    label: "Dry in vacuum",
    icon: FaWind,
    color: "purple",
    summary: ({ values, temperature }) =>
      joinParts([
        withUnit(temperature, "°C"),
        withUnit(values.pressure, "mbar"),
        values.duration,
      ]),
  },
  DRY_WITH_MATERIAL: {
    label: "Dry with material",
    icon: FaCube,
    color: "purple",
  },
  SCAVENGING: {
    label: "Scavenging",
    icon: FaMagnet,
    color: "pink",
    summary: ({ values, materials }) =>
      joinParts([materialsText(materials), values.duration]),
  },
  WAIT: {
    label: "Wait",
    icon: FaClock,
    color: "gray",
    summary: ({ values }) => values.duration,
  },
  STIRRING: {
    label: "Stirring",
    icon: FaSync,
    color: "gray",
    summary: ({ values }) =>
      joinParts([values.duration, withUnit(values.rate, "rpm")]),
  },
  PH_ADJUST: {
    label: "pH adjustment",
    icon: FaBalanceScale,
    color: "yellow",
    summary: ({ values, materials }) =>
      joinParts([
        values.target_ph !== undefined &&
          values.target_ph !== "" &&
          `pH ${values.target_ph}`,
        materialsText(materials),
      ]),
  },
  DISSOLUTION: {
    label: "Dissolution",
    icon: FaVial,
    color: "blue",
  },
  FLASH_CHROMATOGRAPHY: {
    label: "Flash chromatography",
    icon: FaLayerGroup,
    color: "green",
    Component: ChromatographyStep,
    summary: ({ values, materials }) =>
      joinParts([materialsText(materials), values.gradient]),
  },
  OTHER_CHROMATOGRAPHY: {
    label: "Other chromatography",
    icon: FaLayerGroup,
    color: "green",
    Component: ChromatographyStep,
    summary: ({ values, materials }) =>
      joinParts([values.method, materialsText(materials)]),
  },
  DISTILLATION: {
    label: "Distillation",
    icon: FaFire,
    color: "red",
    summary: ({ values, temperature }) =>
      joinParts([
        values.method,
        withUnit(temperature, "°C"),
        withUnit(values.pressure, "mbar"),
      ]),
  },
  CRYSTALLISATION: {
    label: "Crystallisation",
    icon: FaGem,
    color: "green",
    summary: ({ values, materials, temperature }) =>
      joinParts([
        materialsText(materials),
        withUnit(temperature, "°C"),
        values.duration,
      ]),
  },
};

const FALLBACK_STEP_TYPE = {
  icon: FaPencilAlt,
  color: "gray",
};

/**
 * Gets how a step type looks, with a fallback for types the front end doesn't
 * know yet.
 * @param {string} stepType The step type.
 * @returns {object} The step type's label, icon, color, Component and summary.
 */
export const getStepType = (stepType) => {
  const entry = STEP_TYPES_MAP[stepType] ?? {
    ...FALLBACK_STEP_TYPE,
    label: stepType,
  };
  return { Component: GenericStep, summary: defaultSummary, ...entry };
};

/*
 * Which step types the "Add step" menu offers first in each card, in the order
 * chemists usually use them. Every other type is listed after these.
 */
export const COMMON_STEP_TYPES = {
  WORKUP: [
    "ADDITION",
    "EXTRACTION",
    "WASH",
    "DRY_WITH_MATERIAL",
    "FILTRATION",
    "CONCENTRATION",
  ],
  PURIFICATION: [
    "FLASH_CHROMATOGRAPHY",
    "CRYSTALLISATION",
    "DISTILLATION",
    "OTHER_CHROMATOGRAPHY",
    "DRY_IN_VACUUM",
  ],
};
