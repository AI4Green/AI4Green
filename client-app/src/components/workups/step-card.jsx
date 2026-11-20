import {
  Badge,
  Box,
  Button,
  ButtonGroup,
  Collapse,
  FormControl,
  FormLabel,
  HStack,
  Icon,
  IconButton,
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
import { getStepType } from "./step-types";

// how long to wait after the last keystroke before saving
const AUTOSAVE_DELAY_MS = 800;

/** Shows whether the step's latest changes are saved. */
const SaveStatus = ({ status, error }) => {
  if (status === "saving") {
    return (
      <Text fontSize="xs" color="gray.500" role="status">
        Saving…
      </Text>
    );
  }
  if (status === "saved") {
    return (
      <Text fontSize="xs" color="green.600" role="status">
        Saved
      </Text>
    );
  }
  if (status === "error") {
    return (
      <Text fontSize="xs" color="red.600" role="alert">
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
        <IconButton
          size="xs"
          variant="ghost"
          colorScheme="red"
          icon={<FaTrash />}
          aria-label={`Delete ${label} step`}
        />
      </PopoverTrigger>
      <PopoverContent width="auto">
        <PopoverArrow />
        <PopoverBody fontSize="sm">
          Delete this step and anything entered in it?
        </PopoverBody>
        <PopoverFooter display="flex" justifyContent="flex-end">
          <ButtonGroup size="xs">
            <Button variant="ghost" onClick={onClose}>
              Cancel
            </Button>
            <Button
              colorScheme="red"
              onClick={handleDelete}
              isLoading={isDeleting}
            >
              Delete
            </Button>
          </ButtonGroup>
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
 * }} props
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
    defaultIsOpen: false,
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
      borderWidth="1px"
      borderLeftWidth="4px"
      borderLeftColor={`${stepType.color}.400`}
      borderRadius="md"
      bg="white"
    >
      <HStack px={3} py={2} spacing={2}>
        <IconButton
          size="xs"
          variant="ghost"
          icon={isExpanded ? <FaChevronDown /> : <FaChevronRight />}
          aria-label={isExpanded ? "Collapse step" : "Expand step"}
          aria-expanded={isExpanded}
          onClick={onToggle}
        />
        <Badge colorScheme={stepType.color} variant="subtle">
          {position}
        </Badge>
        <Icon as={stepType.icon} color={`${stepType.color}.500`} aria-hidden />
        <Text as="span" fontWeight="semibold" fontSize="sm">
          {stepType.label}
        </Text>
        {!isExpanded && summary && (
          <Text as="span" fontSize="sm" color="gray.600" noOfLines={1}>
            {summary}
          </Text>
        )}
        <Spacer />
        <SaveStatus status={status} error={error} />
        {canEdit && (
          <>
            <IconButton
              size="xs"
              variant="ghost"
              icon={<FaArrowUp />}
              aria-label="Move step up"
              isDisabled={isFirst}
              onClick={() => onMove(-1)}
            />
            <IconButton
              size="xs"
              variant="ghost"
              icon={<FaArrowDown />}
              aria-label="Move step down"
              isDisabled={isLast}
              onClick={() => onMove(1)}
            />
            <DeleteStepButton label={stepType.label} onDelete={onDelete} />
          </>
        )}
      </HStack>

      <Collapse in={isExpanded} animateOpacity>
        <Stack px={3} pb={3} spacing={3}>
          <StepFormContext.Provider value={form}>
            <Component />
          </StepFormContext.Provider>
          <FormControl>
            <FormLabel fontSize="sm" mb={1}>
              Notes
            </FormLabel>
            <Textarea
              size="sm"
              rows={2}
              maxLength={5000}
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
