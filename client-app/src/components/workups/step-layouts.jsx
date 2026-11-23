import {
  Flex,
  FormControl,
  FormLabel,
  Input,
  InputGroup,
  InputRightAddon,
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

// gaps between side-by-side inputs (Chakra spacing units)
const COLUMN_GAP = 6;
const ROW_GAP = 3;

/** The step's temperature in °C (can be below zero). */
const TemperatureField = () => {
  const { definition, temperature, setTemperature, isDisabled } = useStepForm();
  const id = `workup-temperature-${useId()}`;
  return (
    <FormControl id={id} width="150px" maxW="100%">
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
 * Settings (and the temperature, if asked for) side by side on one line, each
 * sized to what it holds, wrapping onto the next line when there isn't room.
 * @param {{ keys?: string[], withTemperature?: boolean }} props `keys` picks
 *   and orders the settings; all of them if left out. The temperature comes
 *   first.
 */
const StepSettings = ({ keys, withTemperature = false }) => {
  const { definition, values, setValue, isDisabled } = useStepForm();
  const fields = keys ? pickFields(definition.fields, keys) : definition.fields;
  const showTemperature = withTemperature && Boolean(definition.temperature);
  if (fields.length === 0 && !showTemperature) return null;

  return (
    <Flex wrap="wrap" align="flex-end" columnGap={COLUMN_GAP} rowGap={ROW_GAP}>
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
    </Flex>
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

/** The material added, then its temperature, addition rate and role. */
export const AdditionStep = () => (
  <Stack spacing={3}>
    <MaterialSlots slots={["material"]} />
    <StepSettings keys={["rate", "role"]} withTemperature />
  </Stack>
);

/** The organic and aqueous phases side by side, then the number of washes. */
export const WashStep = () => (
  <Stack spacing={3}>
    <Flex wrap="wrap" columnGap={COLUMN_GAP + 4} rowGap={ROW_GAP}>
      <MaterialSlots slots={["organic_phase", "aqueous_phase"]} />
    </Flex>
    <StepSettings />
  </Stack>
);

// filter media that use a weighed filter aid, so its mass is asked for
const MEDIA_WITH_FILTER_AID = ["Celite pad", "Silica plug"];

/**
 * The filtration type, medium and what was kept, then the materials. The
 * Celite or silica used is only asked for when the medium needs it, or if
 * one has already been entered (so it isn't hidden when the medium changes).
 */
export const FiltrationStep = () => {
  const { values, materials } = useStepForm();
  const showFilterAid =
    MEDIA_WITH_FILTER_AID.includes(values.filter_medium) ||
    materials.some((row) => row.slot === "filter_aid");
  return (
    <Stack spacing={3}>
      <StepSettings keys={["method", "filter_medium", "phase_kept"]} />
      <MaterialSlots
        slots={
          showFilterAid ? ["filter_aid", "rinse_solvent"] : ["rinse_solvent"]
        }
      />
    </Stack>
  );
};

/** The bath temperature, method and pressure on one line. */
export const ConcentrationStep = () => (
  <StepSettings keys={["method", "pressure"]} withTemperature />
);

/**
 * The method first, then the stationary phase and the eluent (each with its
 * percentage of the mixture), then the gradient.
 */
export const ChromatographyStep = () => (
  <Stack spacing={3}>
    <StepSettings keys={["method"]} />
    <MaterialSlots slots={["stationary_phase", "eluent"]} />
    <StepSettings keys={["gradient"]} />
  </Stack>
);
