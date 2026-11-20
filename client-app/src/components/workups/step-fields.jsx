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
    <FormControl id={id}>
      <FormLabel fontSize="sm" mb={1}>
        {field.label}
      </FormLabel>
      {input}
    </FormControl>
  );
};
