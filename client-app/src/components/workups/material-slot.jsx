import {
  Box,
  Flex,
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

// what a mixture's percentages add up to (matches the backend)
const PROPORTION_TOTAL = 100;

/** Rounds a percentage to 2 decimal places, as the backend stores it. */
const roundPercentage = (value) => Math.round(value * 100) / 100;

/** Whether a material row has a percentage entered. */
const hasProportion = (row) =>
  row.proportion !== "" &&
  row.proportion !== null &&
  row.proportion !== undefined;

/** The sum of a slot's entered percentages. */
const proportionTotal = (rows) =>
  roundPercentage(
    rows.reduce((sum, row) => sum + (Number(row.proportion) || 0), 0),
  );

/**
 * Keeps a mixture's percentages adding up to 100 by setting one row to whatever
 * is left: the last row, or the one before it when the last is being edited
 * (so with two eluents, typing 90 in either sets the other to 10). A single
 * material has no percentage.
 * @param {object[]} rows The slot's material rows, in order.
 * @param {string | null} editedKey The row the user just changed, if any.
 * @returns {object[]} The rows, with the balancing row updated.
 */
const balanceProportions = (rows, editedKey) => {
  if (rows.length < 2) return rows.map((row) => ({ ...row, proportion: "" }));
  const lastIndex = rows.length - 1;
  const balanceIndex =
    rows[lastIndex].key === editedKey ? lastIndex - 1 : lastIndex;
  const others = proportionTotal(
    rows.filter((_, index) => index !== balanceIndex),
  );
  const remainder = roundPercentage(PROPORTION_TOTAL - others);
  // leave it to the user if the others already add up to more than 100
  if (remainder < 0) return rows;
  return rows.map((row, index) =>
    index === balanceIndex ? { ...row, proportion: String(remainder) } : row,
  );
};

/** The total of a mixture's percentages, in red when it isn't 100. */
const ProportionTotal = ({ rows }) => {
  if (!rows.some(hasProportion)) return null;
  const total = proportionTotal(rows);
  const isComplete = total === PROPORTION_TOTAL;
  return (
    <Text
      fontSize="xs"
      color={isComplete ? "gray.600" : "red.600"}
      role={isComplete ? undefined : "alert"}
    >
      Total {total}%{!isComplete && ` (must add up to ${PROPORTION_TOTAL}%)`}
    </Text>
  );
};

/**
 * One material: what it is (from the slot's picker), its percentage of the
 * mixture (for slots with proportions), how much was used and what happened to
 * it afterwards.
 * @param {{
 *   slot: object,
 *   row: object,
 *   label: string,
 *   showProportion: boolean,
 *   onChange: (changes: object) => void,
 *   onRemove?: () => void,
 *   isDisabled: boolean,
 * }} props
 */
const MaterialRow = ({
  slot,
  row,
  label,
  showProportion,
  onChange,
  onRemove,
  isDisabled,
}) => {
  const Picker = slot.picker === "solvent" ? SolventPicker : CompoundPicker;
  const amountInput = (
    <Input
      type="number"
      min={0}
      step="any"
      inputMode="decimal"
      aria-label={`${label} amount`}
      placeholder="Amount"
      width="100px"
      value={row.quantity}
      onChange={(e) => onChange({ quantity: e.target.value })}
      isDisabled={isDisabled}
    />
  );

  // one line: a name box sized for a typical name, then the amount, unit and
  // fate, which move under it together if there isn't room
  return (
    <Flex wrap="wrap" align="center" gap={2}>
      <Box flex="0 0 280px" maxW="100%" minW={0}>
        <Picker
          label={label}
          row={row}
          onChange={onChange}
          isDisabled={isDisabled}
        />
      </Box>
      <HStack spacing={2} flexShrink={0}>
        {showProportion && (
          <InputGroup size="sm" width="auto">
            <Input
              type="number"
              min={0}
              max={PROPORTION_TOTAL}
              step="any"
              inputMode="decimal"
              aria-label={`${label} percentage`}
              title="Percentage of the mixture"
              width="64px"
              value={row.proportion}
              onChange={(e) => onChange({ proportion: e.target.value })}
              isDisabled={isDisabled}
            />
            <InputRightAddon>%</InputRightAddon>
          </InputGroup>
        )}
        <InputGroup size="sm" width="auto">
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
    </Flex>
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
  const showProportions = slot.proportions && rows.length > 1;

  /**
   * Rebalances this slot's percentages within all the step's materials.
   * Only does so once the user has started entering percentages.
   */
  const rebalance = (current, editedKey) => {
    const slotRows = current.filter((item) => item.slot === slot.slot);
    if (!slot.proportions || !slotRows.some(hasProportion)) return current;
    const balanced = new Map(
      balanceProportions(slotRows, editedKey).map((item) => [item.key, item]),
    );
    return current.map((item) => balanced.get(item.key) ?? item);
  };

  const updateRow = (row, changes) =>
    setMaterials((current) => {
      const updated = current.some((item) => item.key === row.key)
        ? current.map((item) =>
            item.key === row.key ? { ...item, ...changes } : item,
          )
        : [...current, { ...row, ...changes }];
      return "proportion" in changes ? rebalance(updated, row.key) : updated;
    });

  const removeRow = (row) =>
    setMaterials((current) =>
      rebalance(
        current.filter((item) => item.key !== row.key),
        null,
      ),
    );

  const addRow = () =>
    setMaterials((current) =>
      rebalance(
        [
          ...current,
          // keep the empty row on screen as the first of the two
          ...(saved.length === 0 ? rows : []),
          emptyMaterial(slot),
        ],
        null,
      ),
    );

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
            showProportion={showProportions}
            onChange={(changes) => updateRow(row, changes)}
            onRemove={saved.includes(row) ? () => removeRow(row) : undefined}
            isDisabled={isDisabled}
          />
          {slot.repeated && <RepeatTotal repeats={repeats} row={row} />}
        </Stack>
      ))}
      {showProportions && <ProportionTotal rows={rows} />}
      {slot.multiple && !isDisabled && (
        <LinkButton
          leftIcon={<FaPlus />}
          alignSelf="flex-start"
          px={0}
          minH="32px"
          onClick={addRow}
        >
          Add another {slot.label.toLowerCase()}
        </LinkButton>
      )}
    </Stack>
  );
};
