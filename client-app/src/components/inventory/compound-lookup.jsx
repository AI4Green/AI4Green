import {
  Box,
  Flex,
  FormControl,
  FormErrorMessage,
  FormHelperText,
  FormLabel,
  Input,
  Spinner,
  Text,
  VStack,
} from "@chakra-ui/react";
import { useBackendApi } from "contexts";
import { useField, useFormikContext } from "formik";
import { useRef, useState } from "react";
import {
  controlProps,
  ErrorMessage,
  HazardTag,
  INV,
  InvButton,
  LinkButton,
  Muted,
} from "./ui.jsx";
import { describeError } from "./use-inventory-action.js";
import { CAS_PATTERN } from "./utils.jsx";

/*
 * Compound picker for Add stock, modelled on the reaction table's reagent
 * lookup. Typing shows matching compounds the workgroup already holds straight
 * away; "Search database" looks the name or CAS number up on the server. A CAS
 * number with no match can be added as a novel compound.
 *
 * Sets the form's `compound` field to the chosen compound:
 * `{ name, cas, hazards, workgroupCompoundId }` for one the workgroup holds,
 * or a lookup result (`{ source, id, name, cas, hazards, workbook,
 * workgroupCompoundId }`) for one it doesn't yet.
 */

const MAX_QUERY_LENGTH = 200;
const MAX_LOCAL_MATCHES = 5;

/**
 * Turns a WorkgroupCompound into the shape the picker stores.
 * @param {object} wc A WorkgroupCompound from the inventory payload.
 * @returns {object} The compound, with its workgroup-level defaults.
 */
export const heldCompound = (wc) => ({
  workgroupCompoundId: wc.id,
  name: wc.compound.name,
  cas: wc.compound.cas,
  hazards: wc.compound.hazards,
  preferredSupplier: wc.preferredSupplier,
  thresholdUnit: wc.thresholdUnit,
});

/**
 * One compound with its CAS number, hazards and where it comes from.
 * @param {object} props
 * @param {object} props.compound The compound to show.
 * @param {React.ReactNode} [props.action] Button shown on the right.
 */
const CompoundCard = ({ compound, action }) => {
  const details = [
    compound.cas || "No CAS number",
    compound.workbook && `Novel compound in ${compound.workbook}`,
    compound.workgroupCompoundId && "Already in this workgroup",
  ].filter(Boolean);

  return (
    <Flex
      justify="space-between"
      align="center"
      gap="12px"
      w="100%"
      px="12px"
      py="8px"
      border="1px solid"
      borderColor={INV.line}
      borderRadius="8px"
      bg="white"
    >
      <Box>
        <Text as="b">{compound.name}</Text>
        <Muted>{details.join(" · ")}</Muted>
        {compound.hazards?.length > 0 && (
          <Box mt="4px">
            {compound.hazards.map((h) => (
              <HazardTag key={h}>{h}</HazardTag>
            ))}
          </Box>
        )}
      </Box>
      {action}
    </Flex>
  );
};

/**
 * The result of a database search: a match, suggestions, an offer to add a
 * new compound, or an error.
 * @param {object} props
 * @param {object} props.lookup The lookup state.
 * @param {(compound: object) => void} props.onSelect Choose a compound.
 * @param {(name: string) => void} props.onSearch Search for a suggested name.
 * @param {(cas: string) => void} props.onAddNew Start adding a new compound.
 */
const LookupResult = ({ lookup, onSelect, onSearch, onAddNew }) => {
  switch (lookup.state) {
    case "loading":
      return (
        <Flex align="center" gap="8px">
          <Spinner size="sm" color={INV.accent} />
          <Text>Searching…</Text>
        </Flex>
      );
    case "error":
      return <ErrorMessage mb={0}>{lookup.message}</ErrorMessage>;
    case "found":
      return (
        <CompoundCard
          compound={lookup.compound}
          action={
            <LinkButton onClick={() => onSelect(lookup.compound)}>
              Use this compound
            </LinkButton>
          }
        />
      );
    case "cas_not_found":
      return (
        <Flex
          align="center"
          justify="space-between"
          wrap="wrap"
          gap="8px"
          w="100%"
        >
          <Text>
            No compound with CAS number <Text as="b">{lookup.cas}</Text> was
            found.
          </Text>
          <InvButton onClick={() => onAddNew(lookup.cas)}>
            Add as new compound
          </InvButton>
        </Flex>
      );
    case "not_found":
      return (
        <Box>
          <Text>
            No exact match.{" "}
            {lookup.suggestions.length > 0
              ? "Did you mean one of these?"
              : "Try the CAS number instead, which also lets you add a new compound."}
          </Text>
          {lookup.suggestions.map((name) => (
            <LinkButton
              key={name}
              px="0"
              mr="12px"
              onClick={() => onSearch(name)}
            >
              {name}
            </LinkButton>
          ))}
        </Box>
      );
    default:
      return null;
  }
};

/**
 * The compound field of the Add stock form.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup adding stock.
 * @param {object[]} props.compounds The workgroup's WorkgroupCompounds.
 * @param {(cas: string) => void} props.onAddNew Opens the new compound form.
 */
export const CompoundLookup = ({ workgroupName, compounds, onAddNew }) => {
  const { inventory } = useBackendApi();
  const { values, setFieldValue } = useFormikContext();
  const [field, meta, helpers] = useField("compound");

  const [query, setQuery] = useState("");
  const [lookup, setLookup] = useState({ state: "idle" });
  // ignores a slow response once a newer search has started
  const latestSearch = useRef(0);

  const selected = field.value;
  const error = meta.touched ? meta.error : undefined;

  /**
   * Chooses a compound. A compound the workgroup already holds also fills in
   * its usual unit and preferred supplier, if those are still empty.
   */
  const select = (compound) => {
    helpers.setValue(compound);
    if (!values.unit && compound.thresholdUnit) {
      setFieldValue("unit", compound.thresholdUnit);
    }
    if (!values.supplier && compound.preferredSupplier) {
      setFieldValue("supplier", compound.preferredSupplier);
    }
  };

  const search = async (text) => {
    const q = text.trim();
    if (!q) return;
    const searchId = ++latestSearch.current;
    setQuery(q);
    setLookup({ state: "loading" });
    try {
      const result = await inventory.lookupCompound(workgroupName, q);
      if (searchId !== latestSearch.current) return;
      setLookup(
        result.status === "found"
          ? { state: "found", compound: result.compound }
          : result.status === "cas_not_found"
          ? { state: "cas_not_found", cas: result.cas }
          : { state: "not_found", suggestions: result.suggestions ?? [] },
      );
    } catch (e) {
      if (searchId !== latestSearch.current) return;
      console.error(e);
      setLookup({ state: "error", message: await describeError(e) });
    }
  };

  const handleQueryChange = (e) => {
    setQuery(e.target.value);
    latestSearch.current += 1; // drop any search still in flight
    setLookup({ state: "idle" });
  };

  // Enter searches; without this it would submit the whole Add stock form
  const handleKeyDown = (e) => {
    if (e.key === "Enter") {
      e.preventDefault();
      search(query);
    }
  };

  const clearSelection = () => {
    helpers.setValue(null);
    setLookup({ state: "idle" });
  };

  if (selected) {
    return (
      <FormControl>
        <FormLabel mb="4px" fontSize="13px" fontWeight={500} color={INV.ink}>
          Compound
        </FormLabel>
        <CompoundCard
          compound={selected}
          action={
            <LinkButton
              aria-label={`Change compound from ${selected.name}`}
              onClick={clearSelection}
            >
              Change
            </LinkButton>
          }
        />
      </FormControl>
    );
  }

  const q = query.trim().toLowerCase();
  const localMatches = q
    ? compounds
        .filter(
          (wc) =>
            wc.compound.name.toLowerCase().includes(q) ||
            wc.compound.cas.includes(q),
        )
        .slice(0, MAX_LOCAL_MATCHES)
    : [];

  return (
    <VStack align="stretch" spacing="8px" w="100%">
      <FormControl id="compound-search" isRequired isInvalid={Boolean(error)}>
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
          Compound
        </FormLabel>
        <Flex gap="8px">
          <Input
            {...controlProps}
            px="12px"
            errorBorderColor={INV.danger}
            type="search"
            placeholder="Name or CAS number"
            autoComplete="off"
            maxLength={MAX_QUERY_LENGTH}
            value={query}
            onChange={handleQueryChange}
            onKeyDown={handleKeyDown}
            onBlur={() => helpers.setTouched(true)}
            autoFocus
          />
          <InvButton
            flexShrink={0}
            onClick={() => search(query)}
            isDisabled={!query.trim()}
          >
            Search database
          </InvButton>
        </Flex>
        {!error && (
          <FormHelperText mt="4px" fontSize="12px" color={INV.muted}>
            {CAS_PATTERN.test(query.trim())
              ? "Search the database for this CAS number. If it isn't there, you can add it."
              : "Search by CAS number to add a compound that isn't in the database yet."}
          </FormHelperText>
        )}
        <FormErrorMessage mt="4px" fontSize="12px" color={INV.danger}>
          {error}
        </FormErrorMessage>
      </FormControl>

      {localMatches.length > 0 && (
        <Box as="section" aria-label="In this workgroup">
          <Muted mb="4px">In this workgroup</Muted>
          <VStack spacing="4px">
            {localMatches.map((wc) => (
              <CompoundCard
                key={wc.id}
                compound={heldCompound(wc)}
                action={
                  <LinkButton onClick={() => select(heldCompound(wc))}>
                    Use
                  </LinkButton>
                }
              />
            ))}
          </VStack>
        </Box>
      )}

      <Box aria-live="polite">
        <LookupResult
          lookup={lookup}
          onSelect={select}
          onSearch={search}
          onAddNew={onAddNew}
        />
      </Box>
    </VStack>
  );
};
