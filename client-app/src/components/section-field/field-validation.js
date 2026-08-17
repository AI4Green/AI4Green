import { isFieldTriggered } from "components/section-field";
import { INPUT_TYPES } from "constants";
import { array, date, mixed, number, object, string } from "yup";

/**
 * Create base validator based on the field type
 * @param {*} fieldType - type of the field
 * @returns
 */
const createBaseValidator = (fieldType) => {
  const {
    Text,
    Number,
    DateAndTime,
    Description,
    FormattedTextInput,
    Multiple,
    Radio,
    ImageFile,
  } = INPUT_TYPES;

  switch (fieldType.toUpperCase()) {
    case Text.toUpperCase():
    case Description.toUpperCase():
    case FormattedTextInput.toUpperCase():
      return string().required("This field is required");

    case Number.toUpperCase():
      return number().required("This field is required");

    case DateAndTime.toUpperCase():
      return date()
        .transform((value) => new Date(value))
        .required("This field is required")
        .nullable();

    case Multiple.toUpperCase():
    case Radio.toUpperCase():
      return array()
        .min(1, "Please select at least one option")
        .required("This field is required");

    case ImageFile.toUpperCase():
      return array()
        .of(
          object().shape({
            file: mixed().when("isNew", {
              is: true,
              then: mixed().required("File is required"),
              otherwise: mixed().notRequired(),
            }),
            caption: string().when("isMarkedForDeletion", {
              is: true,
              then: string().notRequired(),
              otherwise: string().required("Caption is required"),
            }),
          }),
        )
        .min(1, "Please upload at least one image")
        .test(
          "all-marked-for-deletion",
          "Please upload at least one image.",
          (items) => {
            return !items.every((item) => item.isMarkedForDeletion);
          },
        );

    default:
      return mixed().notRequired(); // no validation required for other field types
  }
};

/**
 * Create validation schema for top level fields
 * @param {*} field - field to create schema for
 * @returns
 */
const createTopLevelFieldSchema = (field) => {
  const baseValidator = createBaseValidator(field.inputType.name);
  return [field.id, baseValidator];
};

/**
 * Create validation schema for nested fields
 * @param {*} field - field to create schema for
 * @param {*} allFields - fields collection to search for parent
 * @returns
 */
const createNestedOrTriggeredFieldSchema = (field, allFields) => {
  const parentField = allFields.find((x) => x.triggerField?.id === field.id);
  if (!parentField) return null; // return null if the field has no parent

  let baseValidator = mixed();

  const validator = baseValidator.when(`"${parentField?.id}"`, {
    is: (parentValue) =>
      isFieldTriggered(
        parentField.inputType.name,
        parentField.triggerField.value,
        parentValue,
      ) ?? false, // check if the field is triggered by their parent field
    then: () => {
      baseValidator = createBaseValidator(field.inputType.name);
      return baseValidator.required("This field is required"); // return required validator if the field is triggered
    },
    otherwise: baseValidator.notRequired(),
  });

  return [field.id, validator];
};

/**
 * Create validation schema for the fields collection
 * @param {*} fields - fields collection
 * @returns
 */
export const validationSchema = (fields) => {
  const topLevelFields = fields.filter(
    (field) => !field.hidden && field.mandatory, // assuming hidden fields as not top level
  );
  const nestedOrTriggeredFields = fields.filter(
    (field) => field.hidden && field.mandatory,
  );

  const schemaFields = [
    ...topLevelFields.map((field) => createTopLevelFieldSchema(field)),
    ...nestedOrTriggeredFields
      .map((field) => createNestedOrTriggeredFieldSchema(field, fields))
      .filter(Boolean),
  ];

  return object().shape(Object.fromEntries(schemaFields));
};
