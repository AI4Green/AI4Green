import {
  HStack,
  Input,
  InputGroup,
  InputRightAddon,
  Select,
  Stack,
  Text,
} from "@chakra-ui/react";
import { LinkButton } from "components/inventory/ui";
import { FaPlus, FaTimes } from "react-icons/fa";

import { CompoundPicker, SolventPicker } from "./material-pickers";
import { emptyMaterial, useStepForm } from "./step-form";
import { StepIconButton } from "./step-icon-button";

// matches enums in db
const FATE_LABELS = { WASTE: "Waste", RECOVERED: "Recovered", KEPT: "Kept" };

/** Formats a number without trailing zeros, e.g. 60 rather than 60.0000. */
const formatAmount = (value) =>
  Number(value).toLocaleString(undefined, {
    maximumFractionDigits: 4,
  });

/** "3 × 20 mL = 60 mL", for a slot used once per wash or extraction. */
const RepeatTotal = ({ repeats, row }) => {
  const quantity = Number(row.quantity);
  if (!(repeats > 1) || row.quantity === "" || !Number.isFinite(quantity)) {
    return null;
  }
  return (
    <Text fontSize="xs" color="gray.600">
      {repeats} × {formatAmount(quantity)} {row.unit} ={" "}
      {formatAmount(quantity * repeats)} {row.unit} in total
    </Text>
  );
};

/**
 * One material: what it is (from the slot's picker), how much was used and what
 * happened to it afterwards.
 * @param {{
 *   slot: object,
 *   row: object,
 *   label: string,
 *   onChange: (changes: object) => void,
 *   onRemove?: () => void,
 *   isDisabled: boolean,
 * }} props
 */
const MaterialRow = ({ slot, row, label, onChange, onRemove, isDisabled }) => {
  const Picker = slot.picker === "solvent" ? SolventPicker : CompoundPicker;
  const amountInput = (
    <Input
      type="number"
      min={0}
      step="any"
      inputMode="decimal"
      aria-label={`${label} amount`}
      placeholder="Amount"
      value={row.quantity}
      onChange={(e) => onChange({ quantity: e.target.value })}
      isDisabled={isDisabled}
    />
  );

  return (
    <Stack spacing={1}>
      <Picker
        label={label}
        row={row}
        onChange={onChange}
        isDisabled={isDisabled}
      />
      <HStack spacing={2}>
        <InputGroup size="sm" flex="1">
          {amountInput}
          {slot.units.length === 1 ? (
            <InputRightAddon>{slot.units[0]}</InputRightAddon>
          ) : (
            <Select
              size="sm"
              width="auto"
              aria-label={`${label} unit`}
              value={row.unit}
              onChange={(e) => onChange({ unit: e.target.value })}
              isDisabled={isDisabled}
            >
              {slot.units.map((unit) => (
                <option key={unit} value={unit}>
                  {unit}
                </option>
              ))}
            </Select>
          )}
        </InputGroup>
        <Select
          size="sm"
          width="auto"
          aria-label={`${label} fate`}
          title="What happened to it afterwards"
          value={row.fate}
          onChange={(e) => onChange({ fate: e.target.value })}
          isDisabled={isDisabled}
        >
          {Object.entries(FATE_LABELS).map(([fate, fateLabel]) => (
            <option key={fate} value={fate}>
              {fateLabel}
            </option>
          ))}
        </Select>
        {onRemove && (
          <StepIconButton
            danger
            icon={<FaTimes />}
            aria-label={`Remove ${label}`}
            onClick={onRemove}
            isDisabled={isDisabled}
          />
        )}
      </HStack>
    </Stack>
  );
};

/**
 * A material slot of the step being edited, e.g. a wash's aqueous phase. Shows
 * an empty row until something is entered; slots that take several materials
 * (e.g. an eluent mixture) can add more.
 * @param {{ slot: object }} props The slot's definition from the backend.
 */
export const MaterialSlot = ({ slot }) => {
  const { values, materials, setMaterials, isDisabled } = useStepForm();
  const saved = materials.filter((row) => row.slot === slot.slot);
  // an empty row to type into, added to the step when the user fills it in
  const rows =
    saved.length > 0
      ? saved
      : [{ ...emptyMaterial(slot), key: `blank-${slot.slot}` }];
  const repeats = slot.repeated ? Number(values.repeats) : 0;

  const updateRow = (row, changes) =>
    setMaterials((current) =>
      current.some((item) => item.key === row.key)
        ? current.map((item) =>
            item.key === row.key ? { ...item, ...changes } : item,
          )
        : [...current, { ...row, ...changes }],
    );

  const removeRow = (row) =>
    setMaterials((current) => current.filter((item) => item.key !== row.key));

  return (
    <Stack spacing={2}>
      <Text fontSize="sm" fontWeight="medium">
        {slot.label}
      </Text>
      {rows.map((row, index) => (
        <Stack key={row.key} spacing={1}>
          <MaterialRow
            slot={slot}
            row={row}
            label={index === 0 ? slot.label : `${slot.label} ${index + 1}`}
            onChange={(changes) => updateRow(row, changes)}
            onRemove={saved.includes(row) ? () => removeRow(row) : undefined}
            isDisabled={isDisabled}
          />
          {slot.repeated && <RepeatTotal repeats={repeats} row={row} />}
        </Stack>
      ))}
      {slot.multiple && !isDisabled && (
        <LinkButton
          leftIcon={<FaPlus />}
          alignSelf="flex-start"
          px={0}
          minH="32px"
          onClick={() =>
            setMaterials((current) => [
              ...current,
              // keep the empty row on screen as the first of the two
              ...(saved.length === 0 ? rows : []),
              emptyMaterial(slot),
            ])
          }
        >
          Add another {slot.label.toLowerCase()}
        </LinkButton>
      )}
    </Stack>
  );
};
