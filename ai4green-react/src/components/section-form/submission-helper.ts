import { isFieldTriggered } from "components/section-field";
import { INPUT_TYPES } from "constants";

/**
 * Evaluate field and return the result
 * @param {*} field - field to evaluate
 * @param {*} fields - fields collection
 * @param {*} values - field values collection
 * @returns - returns the evaluated field result
 */
const evaluateField = (field, fields, values) => {
  if (!values[field.id]) return;

  let result = [
    {
      inputType: field.inputType,
      fieldId: field.id,
      fieldResponseId: field.response?.id,
      value: values[field.id],
    },
  ];

  if (!field.triggerField) return result;

  const triggeredFields = fields.filter(
    (childField) => childField.id === field.triggerField.id,
  );

  triggeredFields.forEach((childField) => {
    const triggered = isFieldTriggered(
      field.inputType.name.toUpperCase(),
      field.triggerField.value,
      values[field.id],
    );

    if (triggered) {
      const childResult = evaluateField(childField, fields, values);
      if (childResult) {
        // Concatenate childResult into the result array
        result.push(
          ...(Array.isArray(childResult) ? childResult : [childResult]),
        );
      }
    }
  });

  return result;
};

/**
 * Helper function to prepare the section form data for submission
 * @param {*} fields - fields collection for given section
 * @param {*} values - field values collection
 * @returns - returns the data that is ready for submission
 */
export const prepareSubmissionData = (fields, values) => {
  const data = fields.map((field) => evaluateField(field, fields, values));

  const uniqueData = Object.values(
    data.reduce((uniqueItems, itemArray) => {
      itemArray?.forEach((item) => {
        item && (uniqueItems[item.fieldId] = item);
      });
      return uniqueItems;
    }, {}),
  );

  const {
    fieldResponses,
    newFieldResponses,
    files,
    fileFieldResponses,
    newFiles,
    newFileFieldResponses,
  } = uniqueData.reduce(
    (acc, obj) => {
      return processObject(acc, obj, !obj.fieldResponseId);
    },
    {
      fieldResponses: [],
      newFieldResponses: [],
      files: [],
      fileFieldResponses: [],
      newFiles: [],
      newFileFieldResponses: [],
    },
  );

  return {
    fieldResponses,
    newFieldResponses,
    files,
    fileFieldResponses,
    newFiles,
    newFileFieldResponses,
  };
};

/**
 * Helper function to process the field response object
 * @param {*} acc - accumulator to store the processed data
 * @param {*} obj - field response object
 * @param {*} isNew - flag to determine if the response is new
 * @returns accumulator with processed data
 */
const processObject = (acc, obj, isNew) => {
  const { fieldId, fieldResponseId, inputType } = obj;
  const id = isNew ? fieldId : fieldResponseId;
  const fieldResponseWithFiles = Array.isArray(obj.value)
    ? obj.value
    : obj.value
    ? [obj.value]
    : [];

  if (isFieldFileType(inputType)) {
    processFieldResponseWithFiles(
      acc,
      fieldResponseWithFiles,
      inputType,
      id,
      isNew,
    );
  } else {
    processFieldResponse(acc, obj, id, isNew);
  }

  return acc;
};

/**
 * Helper function to determine if the field's response may include files
 * @param {*} fieldType - field input type
 * @returns true if the field is a file type
 */
const isFieldFileType = (inputType) => {
  const type = inputType.name.toUpperCase();
  return (
    type === INPUT_TYPES.File.toUpperCase() ||
    type === INPUT_TYPES.ImageFile.toUpperCase() ||
    type === INPUT_TYPES.ReactionScheme.toUpperCase()
  );
};

/**
 * Helper function to process field response with files.
 * Ensures files and their field responses are arranged correctly, enabling the backend to relate the files to the correct field response.
 * @param {*} acc - accumulator to store the processed data
 * @param {*} fieldResponseWithFiles - field response with files (array)
 * @param {*} fieldType - field input type
 * @param {*} id - field id or field response id
 * @param {*} isNew - flag to determine if the response is new
 */
const processFieldResponseWithFiles = (
  acc,
  fieldResponseWithFiles,
  inputType,
  id,
  isNew,
) => {
  acc[isNew ? SUBMISSION_KEYS.NewFiles : SUBMISSION_KEYS.Files].push(
    ...fieldResponseWithFiles.map((response) => {
      if (
        inputType.name.toUpperCase() ===
        INPUT_TYPES.ReactionScheme.toUpperCase()
      ) {
        return response.reactionSketch?.reactionImage?.image ?? new Blob();
      }
      return response.isNew ? response.file : new Blob();
    }),
  );

  acc[
    isNew
      ? SUBMISSION_KEYS.NewFieldResponsesWithFiles
      : SUBMISSION_KEYS.FieldResponsesWithFiles
  ].push(
    ...fieldResponseWithFiles.map((response) => ({
      id,
      value: Object.fromEntries(
        Object.entries(response).filter(
          ([key]) => key !== "file" && key !== "image",
        ),
      ),
    })),
  );
};

/**
 * Helper function to process field response
 * @param {*} acc - accumulator to store the processed data
 * @param {*} obj - field response
 * @param {*} id - field id or field response id
 * @param {*} isNew - flag to determine if the response is new
 */
const processFieldResponse = (acc, obj, id, isNew) => {
  acc[
    isNew ? SUBMISSION_KEYS.NewFieldResponses : SUBMISSION_KEYS.FieldResponses
  ].push({
    id,
    value: obj.value,
  });
};

/**
 * Keys used to store the submission data
 */
const SUBMISSION_KEYS = {
  NewFieldResponses: "newFieldResponses",
  FieldResponses: "fieldResponses",
  NewFieldResponsesWithFiles: "newFileFieldResponses",
  NewFiles: "newFiles",
  FieldResponsesWithFiles: "fileFieldResponses",
  Files: "files",
};
