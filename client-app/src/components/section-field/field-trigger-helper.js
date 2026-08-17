import { INPUT_TYPES } from "constants";

/**
 * Helper function to determine whether the field has triggered the child field
 * @param {*} field - field to check
 * @param {*} parentFieldValue - value of the field to check against field trigger value
 * @returns
 */
export const isFieldTriggered = (
  fieldType,
  fieldTriggerValue,
  parentFieldValue,
) => {
  // helper function determining if the trigger condition is met

  const { Text, Description, Multiple, Radio } = INPUT_TYPES;
  switch (fieldType.toUpperCase()) {
    case Text.toUpperCase():
    case Description.toUpperCase():
      return fieldTriggerValue === parentFieldValue;

    case Multiple.toUpperCase():
    case Radio.toUpperCase():
      return parentFieldValue.some((value) => fieldTriggerValue === value.name);

    default:
      return false;
  }
};
