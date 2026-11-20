import {
  Box,
  Collapse,
  FormControl,
  FormLabel,
  HStack,
  Icon,
  Popover,
  PopoverArrow,
  PopoverBody,
  PopoverContent,
  PopoverFooter,
  PopoverTrigger,
  Spacer,
  Stack,
  Text,
  Textarea,
  useDisclosure,
} from "@chakra-ui/react";
import { Badge, focusRing, INV, InvButton } from "components/inventory/ui";
import { useBackendApi } from "contexts";
import { useEffect, useMemo, useRef, useState } from "react";
import {
  FaArrowDown,
  FaArrowUp,
  FaChevronDown,
  FaChevronRight,
  FaTrash,
} from "react-icons/fa";

import { describeWorkupError } from "./errors";
import { StepFormContext, materialPayload, materialRow } from "./step-form";
import { StepIconButton } from "./step-icon-button";
import { getStepType } from "./step-types";

// how long to wait after the last keystroke before saving
const AUTOSAVE_DELAY_MS = 800;

/** Shows whether the step's latest changes are saved. */
const SaveStatus = ({ status, error }) => {
  if (status === "saving") {
    return (
      <Text fontSize="12px" color={INV.muted} role="status">
        Saving…
      </Text>
    );
  }
  if (status === "saved") {
    return (
      <Text fontSize="12px" color={INV.accent} role="status">
        Saved
      </Text>
    );
  }
  if (status === "error") {
    return (
      <Text fontSize="12px" color={INV.danger} role="alert">
        {error}
      </Text>
    );
  }
  return null;
};

/** A delete button that asks for confirmation first. */
const DeleteStepButton = ({ label, onDelete }) => {
  const { isOpen, onOpen, onClose } = useDisclosure();
  const [isDeleting, setIsDeleting] = useState(false);

  const handleDelete = async () => {
    setIsDeleting(true);
    const deleted = await onDelete();
    setIsDeleting(false);
    if (deleted) onClose();
  };

  return (
    <Popover isOpen={isOpen} onOpen={onOpen} onClose={onClose} placement="left">
      <PopoverTrigger>
        <StepIconButton
          danger
          icon={<FaTrash />}
          aria-label={`Delete ${label} step`}
        />
      </PopoverTrigger>
      <PopoverContent
        width="auto"
        borderRadius="12px"
        borderColor={INV.line}
        color={INV.ink}
        _focusVisible={focusRing}
      >
        <PopoverArrow />
        <PopoverBody fontSize="14px">
          Delete this step and anything entered in it?
        </PopoverBody>
        <PopoverFooter
          display="flex"
          justifyContent="flex-end"
          gap="8px"
          borderColor={INV.line}
        >
          <InvButton onClick={onClose} isDisabled={isDeleting}>
            Cancel
          </InvButton>
          <InvButton danger onClick={handleDelete} isLoading={isDeleting}>
            Delete
          </InvButton>
        </PopoverFooter>
      </PopoverContent>
    </Popover>
  );
};

/**
 * One workup step: a header with its type and controls, and the step type's
 * inputs (materials, settings and temperature). Changes are saved
 * automatically shortly after the user stops typing.
 * @param {{
 *   step: object,
 *   position: number,
 *   isFirst: boolean,
 *   isLast: boolean,
 *   definition: { fields: object[], materials: object[], temperature: string | null },
 *   canEdit: boolean,
 *   reactionNames: { workgroup: string, workbook: string, reaction: string },
 *   onMove: (direction: -1 | 1) => Promise<void>,
 *   onDelete: () => Promise<boolean>,
 *   onSaved: (step: object) => void,
 *   defaultExpanded?: boolean,
 * }} props `defaultExpanded` opens the card when it first renders, e.g. for a
 *   step the user has just added; saved steps start collapsed to their summary.
 */
export const StepCard = ({
  step,
  position,
  isFirst,
  isLast,
  definition,
  canEdit,
  reactionNames,
  onMove,
  onDelete,
  onSaved,
  defaultExpanded = false,
}) => {
  const { workups: api } = useBackendApi();
  const stepType = getStepType(step.stepType);
  const { Component } = stepType;

  const [values, setValues] = useState(step.data ?? {});
  const [temperature, setTemperatureState] = useState(step.temperature ?? "");
  const [materials, setMaterialsState] = useState(() =>
    (step.materials ?? []).map(materialRow),
  );
  const [notes, setNotes] = useState(step.notes ?? "");
  const [status, setStatus] = useState("idle");
  const [error, setError] = useState(null);
  const isDirty = useRef(false);
  const { isOpen: isExpanded, onToggle } = useDisclosure({
    defaultIsOpen: defaultExpanded,
  });

  const form = useMemo(
    () => ({
      definition,
      values,
      setValue: (key, value) => {
        isDirty.current = true;
        setValues((current) => ({ ...current, [key]: value }));
      },
      temperature,
      setTemperature: (value) => {
        isDirty.current = true;
        setTemperatureState(value);
      },
      materials,
      setMaterials: (updater) => {
        isDirty.current = true;
        setMaterialsState(updater);
      },
      workbookNames: {
        workgroup: reactionNames.workgroup,
        workbook: reactionNames.workbook,
      },
      isDisabled: !canEdit,
    }),
    [definition, values, temperature, materials, reactionNames, canEdit],
  );

  const handleNotesChange = (e) => {
    isDirty.current = true;
    setNotes(e.target.value);
  };

  useEffect(() => {
    if (!isDirty.current || !canEdit) return undefined;

    const timeout = setTimeout(async () => {
      setStatus("saving");
      try {
        const saved = await api.updateStep(reactionNames, step.id, {
          data: values,
          notes,
          temperature: temperature === "" ? null : temperature,
          materials: materials.map(materialPayload),
        });
        isDirty.current = false;
        setStatus("saved");
        setError(null);
        onSaved(saved);
      } catch (e) {
        console.error(e);
        setStatus("error");
        setError(await describeWorkupError(e));
      }
    }, AUTOSAVE_DELAY_MS);

    return () => clearTimeout(timeout);
    // only the user's edits should trigger a save
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [values, temperature, materials, notes]);

  const summary = stepType.summary({ values, materials, temperature });

  return (
    <Box
      as="section"
      aria-label={`Step ${position}: ${stepType.label}`}
      bg="white"
      color={INV.ink}
      border="1px solid"
      borderColor={INV.line}
      borderLeftWidth="4px"
      borderLeftColor={`${stepType.color}.400`}
      borderRadius="12px"
    >
      <HStack pl="4px" pr="8px" py="2px" spacing="8px">
        <StepIconButton
          icon={isExpanded ? <FaChevronDown /> : <FaChevronRight />}
          aria-label={isExpanded ? "Collapse step" : "Expand step"}
          aria-expanded={isExpanded}
          onClick={onToggle}
        />
        <Badge tone="muted">{position}</Badge>
        <Icon as={stepType.icon} color={`${stepType.color}.600`} aria-hidden />
        <Text as="span" fontWeight={600} fontSize="14px" whiteSpace="nowrap">
          {stepType.label}
        </Text>
        {!isExpanded && summary && (
          <Text as="span" fontSize="13px" color={INV.muted} noOfLines={1}>
            {summary}
          </Text>
        )}
        <Spacer />
        <SaveStatus status={status} error={error} />
        {canEdit && (
          <HStack spacing="0">
            <StepIconButton
              icon={<FaArrowUp />}
              aria-label="Move step up"
              isDisabled={isFirst}
              onClick={() => onMove(-1)}
            />
            <StepIconButton
              icon={<FaArrowDown />}
              aria-label="Move step down"
              isDisabled={isLast}
              onClick={() => onMove(1)}
            />
            <DeleteStepButton label={stepType.label} onDelete={onDelete} />
          </HStack>
        )}
      </HStack>

      <Collapse in={isExpanded} animateOpacity>
        <Stack
          px="16px"
          pt="12px"
          pb="16px"
          spacing="12px"
          borderTop="1px solid"
          borderColor={INV.rowLine}
        >
          <StepFormContext.Provider value={form}>
            <Component />
          </StepFormContext.Provider>
          <FormControl>
            <FormLabel fontSize="13px" fontWeight={500} mb="4px">
              Notes
            </FormLabel>
            <Textarea
              size="sm"
              rows={2}
              maxLength={5000}
              borderRadius="8px"
              borderColor={INV.lineStrong}
              focusBorderColor={INV.accent}
              value={notes}
              onChange={handleNotesChange}
              isDisabled={!canEdit}
            />
          </FormControl>
        </Stack>
      </Collapse>
    </Box>
  );
};
