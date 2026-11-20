import {
  FormControl,
  FormLabel,
  Input,
  InputGroup,
  InputRightAddon,
  SimpleGrid,
  Stack,
} from "@chakra-ui/react";
import { useId } from "react";

import { MaterialSlot } from "./material-slot";
import { StepFieldInput } from "./step-fields";
import { useStepForm } from "./step-form";
import { pickFields } from "./step-utils";

/*
 * Input layouts for step types. A layout says which material slots and settings
 * go where; the inputs read and update the step through useStepForm.
 */

/** The step's temperature in °C (can be below zero). */
const TemperatureField = () => {
  const { definition, temperature, setTemperature, isDisabled } = useStepForm();
  const id = `workup-temperature-${useId()}`;
  return (
    <FormControl id={id}>
      <FormLabel fontSize="sm" mb={1}>
        {definition.temperature}
      </FormLabel>
      <InputGroup size="sm">
        <Input
          type="number"
          step="any"
          inputMode="decimal"
          min={-273.15}
          max={1000}
          value={temperature}
          onChange={(e) => setTemperature(e.target.value)}
          isDisabled={isDisabled}
        />
        <InputRightAddon>°C</InputRightAddon>
      </InputGroup>
    </FormControl>
  );
};

/**
 * Settings (and the temperature, if asked for) in a responsive grid.
 * @param {{ keys?: string[], withTemperature?: boolean }} props `keys` picks
 *   and orders the settings; all of them if left out.
 */
const StepSettings = ({ keys, withTemperature = false }) => {
  const { definition, values, setValue, isDisabled } = useStepForm();
  const fields = keys ? pickFields(definition.fields, keys) : definition.fields;
  const showTemperature = withTemperature && Boolean(definition.temperature);
  if (fields.length === 0 && !showTemperature) return null;

  return (
    <SimpleGrid columns={{ base: 1, md: 2 }} spacing={3}>
      {showTemperature && <TemperatureField />}
      {fields.map((field) => (
        <StepFieldInput
          key={field.key}
          field={field}
          value={values[field.key]}
          onChange={setValue}
          isDisabled={isDisabled}
        />
      ))}
    </SimpleGrid>
  );
};

/**
 * Material slots, one under another.
 * @param {{ slots?: string[] }} props `slots` picks and orders the slots; all
 *   of them if left out.
 */
const MaterialSlots = ({ slots }) => {
  const { definition } = useStepForm();
  const picked = slots
    ? slots
        .map((key) => definition.materials.find((slot) => slot.slot === key))
        .filter(Boolean)
    : definition.materials;
  return picked.map((slot) => <MaterialSlot key={slot.slot} slot={slot} />);
};

/** What was used, then the step's settings and temperature. */
export const GenericStep = () => (
  <Stack spacing={3}>
    <MaterialSlots />
    <StepSettings withTemperature />
  </Stack>
);

/** The method first, then the stationary phase and eluent, then the gradient. */
export const ChromatographyStep = () => (
  <Stack spacing={3}>
    <StepSettings keys={["method"]} />
    <MaterialSlots slots={["stationary_phase", "eluent"]} />
    <StepSettings keys={["gradient"]} />
  </Stack>
);
