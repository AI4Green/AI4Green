import { createContext, useContext } from "react";

/*
 * The state of the step being edited, shared with its layout and inputs so the
 * layouts can just say which fields and material slots go where.
 *
 * Value: {
 *   definition: { fields, materials, temperature },  // from the backend
 *   values, setValue(key, value),                    // settings (WorkupStep.data)
 *   temperature, setTemperature(value),              // °C, as typed
 *   materials, setMaterials(updater),                // material rows
 *   workbookNames: { workgroup, workbook },          // for the lookups
 *   isDisabled,
 * }
 */
export const StepFormContext = createContext(null);

/**
 * Gets the state of the step being edited.
 * @returns {object} See StepFormContext.
 */
export const useStepForm = () => {
  const form = useContext(StepFormContext);
  if (!form) {
    throw new Error("useStepForm must be used inside a StepFormContext");
  }
  return form;
};

let nextMaterialKey = 0;

/**
 * A new, empty material row for a slot. `key` only identifies the row in the
 * browser; it isn't sent to the backend.
 * @param {object} slot The slot definition.
 * @returns {object} The row.
 */
export const emptyMaterial = (slot) => {
  nextMaterialKey += 1;
  return {
    key: `new-${nextMaterialKey}`,
    slot: slot.slot,
    compoundId: null,
    novelCompoundId: null,
    solventName: null,
    name: "",
    quantity: "",
    unit: slot.units[0],
    proportion: "",
    fate: slot.defaultFate,
    flag: null,
  };
};

/**
 * Turns a saved material into a row for editing.
 * @param {object} material A material from the backend.
 * @returns {object} The row.
 */
export const materialRow = (material) => ({
  key: `saved-${material.id}`,
  slot: material.slot,
  compoundId: material.compoundId,
  novelCompoundId: material.novelCompoundId,
  solventName: material.solventName,
  name: material.name,
  quantity: material.quantity ?? "",
  unit: material.unit,
  proportion: material.proportion ?? "",
  fate: material.fate,
  flag: material.flag,
});

/**
 * The material fields the backend accepts (drops the row's browser-only keys).
 * @param {object} row A material row.
 * @returns {object} The material to send.
 */
export const materialPayload = ({
  slot,
  compoundId,
  novelCompoundId,
  solventName,
  name,
  quantity,
  unit,
  fate,
  proportion,
}) => ({
  slot,
  compoundId,
  novelCompoundId,
  solventName,
  name,
  quantity,
  unit,
  fate,
  proportion,
});

/** The link fields of a material that isn't linked to anything. */
export const UNLINKED = {
  compoundId: null,
  novelCompoundId: null,
  solventName: null,
  flag: null,
};
