import {
  FormControl,
  FormErrorMessage,
  FormHelperText,
  FormLabel,
  Input,
  Select,
  Textarea,
} from "@chakra-ui/react";
import { useField } from "formik";
import { controlProps, INV } from "./ui.jsx";

/*
 * Formik-bound form fields in the inventory page's style. Each one shows its
 * validation message once the field has been touched or the form submitted.
 */

/**
 * Label, help text and error message around a single control.
 * @param {object} props
 * @param {string} props.name Formik field name, also used as the control id.
 * @param {string} props.label Label text.
 * @param {boolean} [props.isRequired] Mark the field as required.
 * @param {string} [props.helpText] Hint shown under the control.
 * @param {string} [props.error] Validation message, if invalid.
 */
const FieldFrame = ({ name, label, isRequired, helpText, error, children }) => (
  <FormControl id={name} isRequired={isRequired} isInvalid={Boolean(error)}>
    <FormLabel
      mb="4px"
      fontSize="13px"
      fontWeight={500}
      color={INV.ink}
      requiredIndicator={
        <span aria-hidden="true" style={{ color: INV.danger }}>
          {" "}
          *
        </span>
      }
    >
      {label}
    </FormLabel>
    {children}
    {helpText && !error && (
      <FormHelperText mt="4px" fontSize="12px" color={INV.muted}>
        {helpText}
      </FormHelperText>
    )}
    <FormErrorMessage mt="4px" fontSize="12px" color={INV.danger}>
      {error}
    </FormErrorMessage>
  </FormControl>
);

/**
 * Reads a field's props and its error, shown only after it has been touched.
 * @param {string} name Formik field name.
 * @returns {[object, string | undefined]} The field props and visible error.
 */
const useFieldWithError = (name) => {
  const [field, meta] = useField(name);
  return [field, meta.touched ? meta.error : undefined];
};

/**
 * Single-line text, number or date input.
 * @param {object} props
 * @param {string} props.name Formik field name.
 * @param {string} props.label Label text.
 * @param {boolean} [props.isRequired] Mark the field as required.
 * @param {string} [props.helpText] Hint shown under the input.
 */
export const TextField = ({ name, label, isRequired, helpText, ...props }) => {
  const [field, error] = useFieldWithError(name);
  return (
    <FieldFrame
      name={name}
      label={label}
      isRequired={isRequired}
      helpText={helpText}
      error={error}
    >
      <Input
        {...controlProps}
        px="12px"
        errorBorderColor={INV.danger}
        {...field}
        {...props}
      />
    </FieldFrame>
  );
};

/**
 * Dropdown. Pass the `<option>` elements as children.
 * @param {object} props
 * @param {string} props.name Formik field name.
 * @param {string} props.label Label text.
 * @param {boolean} [props.isRequired] Mark the field as required.
 * @param {string} [props.helpText] Hint shown under the dropdown.
 */
export const SelectField = ({
  name,
  label,
  isRequired,
  helpText,
  children,
  ...props
}) => {
  const [field, error] = useFieldWithError(name);
  return (
    <FieldFrame
      name={name}
      label={label}
      isRequired={isRequired}
      helpText={helpText}
      error={error}
    >
      <Select
        {...controlProps}
        ps="12px"
        errorBorderColor={INV.danger}
        {...field}
        {...props}
      >
        {children}
      </Select>
    </FieldFrame>
  );
};

/**
 * Multi-line text input.
 * @param {object} props
 * @param {string} props.name Formik field name.
 * @param {string} props.label Label text.
 * @param {boolean} [props.isRequired] Mark the field as required.
 * @param {string} [props.helpText] Hint shown under the textarea.
 */
export const TextAreaField = ({
  name,
  label,
  isRequired,
  helpText,
  ...props
}) => {
  const [field, error] = useFieldWithError(name);
  return (
    <FieldFrame
      name={name}
      label={label}
      isRequired={isRequired}
      helpText={helpText}
      error={error}
    >
      <Textarea
        {...controlProps}
        h="auto"
        minH="88px"
        px="12px"
        py="10px"
        errorBorderColor={INV.danger}
        {...field}
        {...props}
      />
    </FieldFrame>
  );
};
