import { Badge, Box, HStack, Input, Text } from "@chakra-ui/react";
import { useWorkupSolvents } from "api";
import { useBackendApi } from "contexts";
import { useEffect, useId, useRef, useState } from "react";

import { flagLabel, flagStyle } from "./chem21";
import { UNLINKED, useStepForm } from "./step-form";

// most options shown in a picker's list at once
const MAX_OPTIONS = 50;
// how long to wait after the last keystroke before looking a compound up
const LOOKUP_DELAY_MS = 400;
// shortest text worth looking up
const MIN_LOOKUP_LENGTH = 2;

/** Whether a material row is linked to a compound, novel compound or solvent. */
const isLinked = (row) =>
  Boolean(row.compoundId || row.novelCompoundId || row.solventName);

/** The link fields for a compound found by the lookup. */
const compoundLink = (compound) => ({
  ...UNLINKED,
  name: compound.name,
  compoundId: compound.source === "compound" ? compound.id : null,
  novelCompoundId: compound.source === "novel" ? compound.id : null,
});

/**
 * A text input with a dropdown list of options, following the ARIA combobox
 * pattern. The list is rendered in place (not in a portal), so it stays inside
 * the embedded React root on the reaction page.
 * @param {{
 *   label: string,
 *   text: string,
 *   onTextChange: (text: string) => void,
 *   options: { key: string, label: string, detail?: string, style?: object }[],
 *   onSelect: (option: object) => void,
 *   onCommit?: () => void,
 *   isDisabled: boolean,
 * }} props `onCommit` runs when the input loses focus.
 */
const MaterialCombobox = ({
  label,
  text,
  onTextChange,
  options,
  onSelect,
  onCommit,
  isDisabled,
}) => {
  const [isOpen, setIsOpen] = useState(false);
  const [activeIndex, setActiveIndex] = useState(-1);
  const listId = useId();
  const showList = isOpen && options.length > 0;

  const select = (option) => {
    onSelect(option);
    setIsOpen(false);
    setActiveIndex(-1);
  };

  const handleKeyDown = (e) => {
    if (e.key === "ArrowDown") {
      e.preventDefault();
      setIsOpen(true);
      setActiveIndex((index) => Math.min(index + 1, options.length - 1));
    } else if (e.key === "ArrowUp") {
      e.preventDefault();
      setActiveIndex((index) => Math.max(index - 1, 0));
    } else if (e.key === "Enter" && showList && options[activeIndex]) {
      e.preventDefault();
      select(options[activeIndex]);
    } else if (e.key === "Escape") {
      setIsOpen(false);
    }
  };

  return (
    <Box position="relative" flex="1" minW={0}>
      <Input
        size="sm"
        role="combobox"
        aria-label={label}
        aria-autocomplete="list"
        aria-expanded={showList}
        aria-controls={listId}
        aria-activedescendant={
          showList && activeIndex >= 0 ? `${listId}-${activeIndex}` : undefined
        }
        autoComplete="off"
        maxLength={200}
        value={text}
        onChange={(e) => {
          onTextChange(e.target.value);
          setIsOpen(true);
          setActiveIndex(-1);
        }}
        onFocus={() => setIsOpen(true)}
        onBlur={() => {
          setIsOpen(false);
          onCommit?.();
        }}
        onKeyDown={handleKeyDown}
        isDisabled={isDisabled}
      />
      {showList && (
        <Box
          as="ul"
          role="listbox"
          id={listId}
          aria-label={`${label} options`}
          position="absolute"
          zIndex="dropdown"
          top="100%"
          left={0}
          right={0}
          mt={1}
          maxH="240px"
          overflowY="auto"
          bg="white"
          borderWidth="1px"
          borderRadius="md"
          boxShadow="md"
          listStyleType="none"
          m={0}
          p={1}
        >
          {options.map((option, index) => (
            <Box
              as="li"
              key={option.key}
              id={`${listId}-${index}`}
              role="option"
              aria-selected={index === activeIndex}
              // keep focus in the input, so choosing doesn't count as leaving it
              onMouseDown={(e) => {
                e.preventDefault();
                select(option);
              }}
              cursor="pointer"
              px={2}
              py={1}
              mb="1px"
              borderRadius="sm"
              fontSize="sm"
              outline={index === activeIndex ? "2px solid" : "none"}
              outlineColor="blue.400"
              display="flex"
              justifyContent="space-between"
              gap={2}
              {...(option.style ?? { _hover: { bg: "gray.100" } })}
            >
              <span>{option.label}</span>
              {option.detail && (
                <Text as="span" fontSize="xs" opacity={0.8}>
                  {option.detail}
                </Text>
              )}
            </Box>
          ))}
        </Box>
      )}
    </Box>
  );
};

/** Shows what a material is linked to: its CHEM21 flag, or whether it was found. */
const LinkIndicator = ({ row, colours }) => {
  if (row.flag) {
    return (
      <Badge {...flagStyle(row.flag, colours)} borderWidth="1px" fontSize="2xs">
        {flagLabel(row.flag)}
      </Badge>
    );
  }
  if (isLinked(row)) {
    return (
      <Badge colorScheme="green" variant="subtle" fontSize="2xs">
        In database
      </Badge>
    );
  }
  if (row.name.trim()) {
    return (
      <Text fontSize="2xs" color="gray.500" whiteSpace="nowrap">
        Not in database
      </Text>
    );
  }
  return null;
};

/**
 * Looks typed compound names up after a pause, ignoring results for text the
 * user has since changed.
 * @param {(result: object, text: string) => void} onResult Called with each result.
 * @returns {{ schedule: (text: string) => void, lookupNow: (text: string) => void }}
 */
const useCompoundLookup = (onResult) => {
  const { workbookNames } = useStepForm();
  const { workups: api } = useBackendApi();
  const latestText = useRef("");
  const timer = useRef(null);

  useEffect(() => () => clearTimeout(timer.current), []);

  const lookupNow = async (text) => {
    clearTimeout(timer.current);
    latestText.current = text;
    if (text.length < MIN_LOOKUP_LENGTH) return;
    try {
      const result = await api.lookupCompound(workbookNames, text);
      if (latestText.current === text) onResult(result, text);
    } catch (e) {
      // a failed lookup just leaves the material unlinked
      console.error(e);
    }
  };

  const schedule = (text) => {
    clearTimeout(timer.current);
    latestText.current = text;
    if (text.length < MIN_LOOKUP_LENGTH) return;
    timer.current = setTimeout(() => lookupNow(text), LOOKUP_DELAY_MS);
  };

  return { schedule, lookupNow };
};

/**
 * Picks a solvent from the workbook's CHEM21 solvent list, coloured by each
 * solvent's CHEM21 flag in the user's colours, as on the reaction table. Text
 * that isn't on the list is looked up in the compound database when the input
 * loses focus, and kept as an unlinked name if it isn't found there either.
 * @param {{ label: string, row: object, onChange: (changes: object) => void, isDisabled: boolean }} props
 */
export const SolventPicker = ({ label, row, onChange, isDisabled }) => {
  const { workbookNames } = useStepForm();
  const { data } = useWorkupSolvents(workbookNames);
  const solvents = data?.solvents ?? [];
  const colours = data?.colours;
  const { lookupNow } = useCompoundLookup((result) => {
    if (result.status === "found") onChange(compoundLink(result.compound));
  });

  const query = row.name.trim().toLowerCase();
  const matches = (
    isLinked(row) || !query
      ? solvents
      : solvents.filter((solvent) => solvent.name.toLowerCase().includes(query))
  ).slice(0, MAX_OPTIONS);

  const choose = (solvent) =>
    onChange({
      name: solvent.name,
      compoundId: solvent.compoundId,
      novelCompoundId: solvent.novelCompoundId,
      solventName: solvent.solventName,
      flag: solvent.flag,
    });

  const handleCommit = () => {
    const text = row.name.trim();
    if (!text || isLinked(row)) return;
    const exact = solvents.find(
      (solvent) => solvent.name.toLowerCase() === text.toLowerCase(),
    );
    if (exact) {
      choose(exact);
    } else {
      lookupNow(text);
    }
  };

  return (
    <HStack spacing={2} align="center">
      <MaterialCombobox
        label={label}
        text={row.name}
        onTextChange={(text) => onChange({ ...UNLINKED, name: text })}
        options={matches.map((solvent) => ({
          key: solvent.name,
          label: solvent.name,
          detail: flagLabel(solvent.flag),
          style: flagStyle(solvent.flag, colours),
          solvent,
        }))}
        onSelect={(option) => choose(option.solvent)}
        onCommit={handleCommit}
        isDisabled={isDisabled}
      />
      <LinkIndicator row={row} colours={colours} />
    </HStack>
  );
};

/**
 * Picks a compound by name or CAS number from the compound database and the
 * user's novel compounds (the inventory's lookup). An exact match links the
 * material; otherwise similar names are suggested, and text that matches
 * nothing is kept as an unlinked name.
 * @param {{ label: string, row: object, onChange: (changes: object) => void, isDisabled: boolean }} props
 */
export const CompoundPicker = ({ label, row, onChange, isDisabled }) => {
  const [suggestions, setSuggestions] = useState([]);
  const { schedule, lookupNow } = useCompoundLookup((result) => {
    if (result.status === "found") {
      setSuggestions([]);
      onChange(compoundLink(result.compound));
    } else {
      setSuggestions(result.suggestions ?? []);
    }
  });

  const handleTextChange = (text) => {
    onChange({ ...UNLINKED, name: text });
    setSuggestions([]);
    schedule(text.trim());
  };

  const handleSelect = ({ label: name }) => {
    onChange({ ...UNLINKED, name });
    setSuggestions([]);
    lookupNow(name);
  };

  return (
    <HStack spacing={2} align="center">
      <MaterialCombobox
        label={label}
        text={row.name}
        onTextChange={handleTextChange}
        options={suggestions.map((name) => ({ key: name, label: name }))}
        onSelect={handleSelect}
        isDisabled={isDisabled}
      />
      <LinkIndicator row={row} />
    </HStack>
  );
};
