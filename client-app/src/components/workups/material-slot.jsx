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
import { useState } from "react";
import { FaPlus, FaTimes } from "react-icons/fa";

import { CompoundPicker, SolventPicker } from "./material-pickers";
import {
  PROPORTION_TOTAL,
  balanceProportions,
  hasProportion,
  isProportionComplete,
  mixtureVolume,
  proportionTotal,
  splitVolume,
} from "./mixture";
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
 * An amount with its unit: a fixed unit, or a choice when the slot allows
 * several.
 * @param {{
 *   label: string,
 *   value: string,
 *   unit: string,
 *   units: string[],
 *   onValueChange: (value: string) => void,
 *   onUnitChange: (unit: string) => void,
 *   isDisabled: boolean,
 * }} props `label` names the inputs for screen readers.
 */
const AmountInput = ({
  label,
  value,
  unit,
  units,
  onValueChange,
  onUnitChange,
  isDisabled,
}) => (
  <InputGroup size="sm" width="auto">
    <Input
      type="number"
      min={0}
      step="any"
      inputMode="decimal"
      aria-label={label}
      placeholder="Amount"
      width="100px"
      value={value}
      onChange={(e) => onValueChange(e.target.value)}
      isDisabled={isDisabled}
    />
    {units.length === 1 ? (
      <InputRightAddon>{units[0]}</InputRightAddon>
    ) : (
      <Select
        size="sm"
        width="auto"
        aria-label={`${label} unit`}
        value={unit}
        onChange={(e) => onUnitChange(e.target.value)}
        isDisabled={isDisabled}
      >
        {units.map((option) => (
          <option key={option} value={option}>
            {option}
          </option>
        ))}
      </Select>
    )}
  </InputGroup>
);

/**
 * One material: what it is (from the slot's picker), how much was used and what
 * happened to it afterwards. In a mixture it has a percentage instead of its
 * own amount, which comes from the mixture's total volume.
 * @param {{
 *   slot: object,
 *   row: object,
 *   label: string,
 *   isMixture: boolean,
 *   onChange: (changes: object) => void,
 *   onRemove?: () => void,
 *   isDisabled: boolean,
 * }} props
 */
const MaterialRow = ({
  slot,
  row,
  label,
  isMixture,
  onChange,
  onRemove,
  isDisabled,
}) => {
  const Picker = slot.picker === "solvent" ? SolventPicker : CompoundPicker;

  // one line: a name box sized for a typical name, then the amount (or
  // percentage) and fate, which move under it together if there isn't room
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
        {isMixture ? (
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
        ) : (
          <AmountInput
            label={`${label} amount`}
            value={row.quantity}
            unit={row.unit}
            units={slot.units}
            onValueChange={(quantity) => onChange({ quantity })}
            onUnitChange={(unit) => onChange({ unit })}
            isDisabled={isDisabled}
          />
        )}
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
      Total {total}%
      {!isComplete &&
        ` (must add up to ${PROPORTION_TOTAL}% to split the volume)`}
    </Text>
  );
};

/**
 * The volume of a mixture as made up, e.g. 500 mL of 90:10 hexane/EtOAc,
 * with each material's share of it once the percentages add up to 100.
 * @param {{
 *   slot: object,
 *   rows: object[],
 *   volume: string,
 *   onVolumeChange: (volume: string) => void,
 *   onUnitChange: (unit: string) => void,
 *   isDisabled: boolean,
 * }} props
 */
const MixtureVolume = ({
  slot,
  rows,
  volume,
  onVolumeChange,
  onUnitChange,
  isDisabled,
}) => {
  const unit = rows[0].unit;
  const shares = isProportionComplete(rows)
    ? rows
        .filter((row) => row.name && row.quantity !== "")
        .map((row) => `${row.name} ${formatAmount(row.quantity)} ${unit}`)
    : [];
  return (
    <Stack spacing={1}>
      <HStack spacing={3} flexWrap="wrap" align="center">
        {/* a span, so the reaction page's paragraph margin doesn't push it up */}
        <Text as="span" fontSize="sm" lineHeight="32px" m={0}>
          Total volume used:
        </Text>
        <AmountInput
          label={`${slot.label} total volume`}
          value={volume}
          unit={unit}
          units={slot.units}
          onValueChange={onVolumeChange}
          onUnitChange={onUnitChange}
          isDisabled={isDisabled}
        />
      </HStack>
      <ProportionTotal rows={rows} />
      {shares.length > 0 && (
        <Text fontSize="xs" color="gray.600">
          {shares.join(" + ")}
        </Text>
      )}
    </Stack>
  );
};

/**
 * A material slot of the step being edited, e.g. a wash's aqueous phase. Shows
 * an empty row until something is entered; slots that take several materials
 * (e.g. an eluent mixture) can add more. In a slot with proportions, two or
 * more materials make a mixture: each gets a percentage, and one total volume
 * is split between them.
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
  const isMixture = Boolean(slot.proportions) && rows.length > 1;
  // the mixture's volume as typed; each material's amount is worked out from it
  const [volume, setVolume] = useState(() =>
    isMixture ? mixtureVolume(saved) : "",
  );

  /** Applies a change to this slot's rows within all the step's materials. */
  const updateSlotRows = (current, transform) => {
    const slotRows = current.filter((item) => item.slot === slot.slot);
    const updated = new Map(
      transform(slotRows).map((item) => [item.key, item]),
    );
    return current.map((item) => updated.get(item.key) ?? item);
  };

  /**
   * Rebalances the percentages (once the user has started entering them) and
   * splits the mixture's volume between its materials.
   */
  const settleMixture = (current, editedKey, mixtureVolumeText) => {
    if (!slot.proportions) return current;
    return updateSlotRows(current, (slotRows) => {
      const balanced = slotRows.some(hasProportion)
        ? balanceProportions(slotRows, editedKey)
        : slotRows;
      return balanced.length > 1
        ? splitVolume(balanced, mixtureVolumeText)
        : balanced;
    });
  };

  const updateRow = (row, changes) =>
    setMaterials((current) => {
      const updated = current.some((item) => item.key === row.key)
        ? current.map((item) =>
            item.key === row.key ? { ...item, ...changes } : item,
          )
        : [...current, { ...row, ...changes }];
      return "proportion" in changes
        ? settleMixture(updated, row.key, volume)
        : updated;
    });

  const removeRow = (row) =>
    setMaterials((current) => {
      const remaining = current.filter((item) => item.key !== row.key);
      const settled = settleMixture(remaining, null, volume);
      // the last material left is the whole of what was the mixture
      if (isMixture && saved.length === 2 && volume !== "") {
        return updateSlotRows(settled, (slotRows) =>
          slotRows.map((item) => ({ ...item, quantity: volume })),
        );
      }
      return settled;
    });

  const addRow = () => {
    // a single material's amount becomes the new mixture's volume
    const nextVolume = isMixture ? volume : mixtureVolume(saved);
    setVolume(nextVolume);
    setMaterials((current) =>
      settleMixture(
        [
          ...current,
          // keep the empty row on screen as the first of the two
          ...(saved.length === 0 ? rows : []),
          // in the same unit, as a mixture has one volume
          { ...emptyMaterial(slot), unit: rows[0].unit },
        ],
        null,
        nextVolume,
      ),
    );
  };

  const changeVolume = (text) => {
    setVolume(text);
    setMaterials((current) =>
      updateSlotRows(current, (slotRows) => splitVolume(slotRows, text)),
    );
  };

  const changeMixtureUnit = (unit) =>
    setMaterials((current) =>
      updateSlotRows(current, (slotRows) =>
        slotRows.map((item) => ({ ...item, unit })),
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
            isMixture={isMixture}
            onChange={(changes) => updateRow(row, changes)}
            onRemove={saved.includes(row) ? () => removeRow(row) : undefined}
            isDisabled={isDisabled}
          />
          {slot.repeated && <RepeatTotal repeats={repeats} row={row} />}
        </Stack>
      ))}
      {isMixture && (
        <MixtureVolume
          slot={slot}
          rows={rows}
          volume={volume}
          onVolumeChange={changeVolume}
          onUnitChange={changeMixtureUnit}
          isDisabled={isDisabled}
        />
      )}
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
