import {
  FormControl,
  FormLabel,
  Input,
  InputGroup,
  InputRightAddon,
  Select,
} from "@chakra-ui/react";
import { useId } from "react";

/*
 * Inputs for the settings a step type records. The field definitions come from
 * the backend (`/api/workups/step-definitions`): each has a `key`, `label`,
 * `kind` ("text" | "number" | "choice") and optionally a `unit`, `choices` and
 * `maximum`. Values are kept as the user typed them; the backend cleans and
 * validates them when the step is saved.
 */

/*
 * How wide each kind of field is, so a short value (e.g. the number of washes)
 * gets a small box rather than stretching across the card.
 */
const FIELD_WIDTHS = {
  number: "100px",
  numberWithUnit: "150px",
  choice: "190px",
  text: "240px",
};

/**
 * The width of a field's input.
 * @param {object} field The field definition.
 * @returns {string} A CSS width.
 */
const fieldWidth = (field) => {
  if (field.kind === "number") {
    return field.unit ? FIELD_WIDTHS.numberWithUnit : FIELD_WIDTHS.number;
  }
  return FIELD_WIDTHS[field.kind] ?? FIELD_WIDTHS.text;
};

/**
 * One field's input, chosen by its kind.
 * @param {{ field: object, value: *, onChange: (key: string, value: *) => void, isDisabled: boolean }} props
 */
export const StepFieldInput = ({ field, value, onChange, isDisabled }) => {
  const id = `workup-${field.key}-${useId()}`;
  const current = value ?? "";

  let input;
  if (field.kind === "choice") {
    input = (
      <Select
        size="sm"
        placeholder="Select…"
        value={current}
        onChange={(e) => onChange(field.key, e.target.value)}
        isDisabled={isDisabled}
      >
        {field.choices.map((choice) => (
          <option key={choice} value={choice}>
            {choice}
          </option>
        ))}
      </Select>
    );
  } else if (field.kind === "number") {
    const numberInput = (
      <Input
        size="sm"
        type="number"
        min={0}
        max={field.maximum ?? undefined}
        step="any"
        inputMode="decimal"
        value={current}
        onChange={(e) => onChange(field.key, e.target.value)}
        isDisabled={isDisabled}
      />
    );
    input = field.unit ? (
      <InputGroup size="sm">
        {numberInput}
        <InputRightAddon>{field.unit}</InputRightAddon>
      </InputGroup>
    ) : (
      numberInput
    );
  } else {
    input = (
      <Input
        size="sm"
        maxLength={500}
        value={current}
        onChange={(e) => onChange(field.key, e.target.value)}
        isDisabled={isDisabled}
      />
    );
  }

  return (
    <FormControl id={id} width={fieldWidth(field)} maxW="100%">
      <FormLabel fontSize="sm" mb={1}>
        {field.label}
      </FormLabel>
      {input}
    </FormControl>
  );
};
