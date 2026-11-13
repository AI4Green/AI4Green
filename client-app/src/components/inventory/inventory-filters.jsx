import {
  Button,
  Checkbox,
  Flex,
  FormControl,
  FormLabel,
  Input,
  Select,
} from "@chakra-ui/react";
import { STRUCTURE_MODE_LABELS } from "components/core/structure-sketcher/structure-search-labels.js";
import {
  controlProps,
  focusRing,
  INV,
  InvButton,
  LinkButton,
  underline,
} from "./ui.jsx";

/**
 * Labelled filter control.
 * @param {object} props
 * @param {string} props.id Id for the control, linked to the label.
 * @param {string} props.label Label text.
 * @param {boolean} [props.grow] Take up the remaining row width.
 */
const Field = ({ id, label, grow = false, children }) => (
  <FormControl
    id={id}
    display="flex"
    flexDirection="column"
    gap="4px"
    w="auto"
    flex={grow ? 1 : undefined}
    minW={grow ? "260px" : undefined}
  >
    <FormLabel
      m={0}
      fontSize="12px"
      fontWeight={500}
      color={INV.muted}
      lineHeight="inherit"
    >
      {label}
    </FormLabel>
    {children}
  </FormControl>
);

/**
 * One half of a two-way toggle, such as "By compound / All containers". The
 * selected half is underlined as well as coloured, so the state doesn't rely
 * on colour alone.
 */
export const SegmentButton = ({ isPressed, position, ...props }) => (
  <Button
    variant="unstyled"
    display="inline-flex"
    h="auto"
    minW={0}
    minH="44px"
    px="16px"
    lineHeight="inherit"
    fontSize="14px"
    fontWeight={400}
    bg="white"
    color={INV.ink}
    border="1px solid"
    borderColor={INV.control}
    borderRadius={position === "first" ? "8px 0 0 8px" : "0 8px 8px 0"}
    borderLeftWidth={position === "last" ? 0 : "1px"}
    aria-pressed={isPressed}
    _pressed={{
      bg: INV.accentSoft,
      borderColor: INV.accent,
      color: INV.accent,
      fontWeight: 500,
      ...underline,
    }}
    _hover={underline}
    _focusVisible={{ ...focusRing, ...underline }}
    {...props}
  />
);

/**
 * Shows the structure search in use, with a button to remove it.
 * @param {object} props
 * @param {{ ids: Set<number>, mode: string }} props.structure The search result.
 * @param {() => void} props.onEdit Opens the structure search again.
 * @param {() => void} props.onClear Removes the structure filter.
 */
const StructureFilterChip = ({ structure, onEdit, onClear }) => (
  <Flex
    role="status"
    align="center"
    gap="4px"
    minH="44px"
    ps="12px"
    border="1px solid"
    borderColor={INV.accent}
    borderRadius="8px"
    bg={INV.accentSoft}
    fontSize="14px"
  >
    Structure: {STRUCTURE_MODE_LABELS[structure.mode].toLowerCase()} (
    {structure.ids.size} {structure.ids.size === 1 ? "compound" : "compounds"})
    <LinkButton onClick={onEdit}>Change</LinkButton>
    <LinkButton onClick={onClear} aria-label="Clear the structure search">
      Clear
    </LinkButton>
  </Flex>
);

export default function InventoryFilters({
  filters,
  onChange,
  locations,
  suppliers,
  view,
  onViewChange,
  onOpenStructureSearch,
}) {
  const set = (key) => (e) =>
    onChange({
      ...filters,
      [key]: e.target.type === "checkbox" ? e.target.checked : e.target.value,
    });

  return (
    <Flex wrap="wrap" gap="16px" align="flex-end">
      <Field id="inv-search" label="Search" grow>
        <Input
          {...controlProps}
          type="search"
          px="12px"
          placeholder="Name, CAS number, barcode or lot"
          value={filters.search}
          onChange={set("search")}
        />
      </Field>

      {filters.structure ? (
        <StructureFilterChip
          structure={filters.structure}
          onEdit={onOpenStructureSearch}
          onClear={() => onChange({ ...filters, structure: null })}
        />
      ) : (
        <InvButton onClick={onOpenStructureSearch}>
          Search by structure
        </InvButton>
      )}

      <Field id="inv-location" label="Location">
        <Select
          {...controlProps}
          ps="12px"
          value={filters.location}
          onChange={set("location")}
        >
          <option value="">All locations</option>
          {locations.map((l) => (
            <option key={l} value={l}>
              {l}
            </option>
          ))}
        </Select>
      </Field>

      <Field id="inv-supplier" label="Supplier">
        <Select
          {...controlProps}
          ps="12px"
          value={filters.supplier}
          onChange={set("supplier")}
        >
          <option value="">All suppliers</option>
          {suppliers.map((s) => (
            <option key={s} value={s}>
              {s}
            </option>
          ))}
        </Select>
      </Field>

      <Field id="inv-status" label="Status">
        <Select
          {...controlProps}
          ps="12px"
          value={filters.status}
          onChange={set("status")}
        >
          <option value="">Any status</option>
          <option value="low">Low stock</option>
          <option value="expiring">Expiring soon</option>
        </Select>
      </Field>

      <Flex role="group" aria-label="View">
        <SegmentButton
          position="first"
          isPressed={view === "compound"}
          onClick={() => onViewChange("compound")}
        >
          By compound
        </SegmentButton>
        <SegmentButton
          position="last"
          isPressed={view === "containers"}
          // mockup draws an accent seam when the left half is selected
          borderLeftWidth={view === "compound" ? "1px" : 0}
          borderLeftColor={view === "compound" ? INV.accent : undefined}
          onClick={() => onViewChange("containers")}
        >
          All containers
        </SegmentButton>
      </Flex>

      <Checkbox
        isChecked={filters.showArchived}
        onChange={set("showArchived")}
        spacing="8px"
        minH="44px"
        sx={{
          ".chakra-checkbox__label": { fontSize: "14px" },
          ".chakra-checkbox__control": { borderColor: INV.control },
          ".chakra-checkbox__control[data-checked], .chakra-checkbox__control[data-checked][data-hover]":
            { bg: INV.accent, borderColor: INV.accent },
          ".chakra-checkbox__control[data-focus-visible]": focusRing,
          ".chakra-checkbox__control[data-focus-visible] ~ .chakra-checkbox__label":
            underline,
        }}
      >
        Show archived
      </Checkbox>
    </Flex>
  );
}
