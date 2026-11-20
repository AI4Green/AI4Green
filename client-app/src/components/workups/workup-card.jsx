import {
  Alert,
  AlertIcon,
  Box,
  Button,
  HStack,
  Menu,
  MenuButton,
  MenuDivider,
  MenuGroup,
  MenuItem,
  MenuList,
  Spinner,
  Stack,
  Text,
  useDisclosure,
} from "@chakra-ui/react";
import { useStepDefinitions, useWorkupSteps, useWorkupTemplates } from "api";
import { useBackendApi } from "contexts";
import { useEffect, useMemo, useState } from "react";
import { FaChevronDown, FaCopy, FaPlus, FaSave } from "react-icons/fa";

import { describeWorkupError } from "./errors";
import { StepCard } from "./step-card";
import { StepTypeLabel } from "./step-type-label";
import { COMMON_STEP_TYPES } from "./step-types";
import { SaveTemplateModal, TemplatesModal } from "./template-modals";

const CATEGORY_LABELS = { WORKUP: "workup", PURIFICATION: "purification" };
// how long a success message stays in the card
const NOTICE_DURATION_MS = 4000;
// for a step type the backend no longer defines: only its notes can be edited
const EMPTY_DEFINITION = { fields: [], materials: [], temperature: null };

/** The "Add step" menu: this card's common step types first, then the rest. */
const AddStepMenu = ({ category, stepTypes, onAdd, isLoading }) => {
  const common = COMMON_STEP_TYPES[category] ?? [];
  const others = stepTypes.filter(
    (type) => !common.includes(type) && type !== "CUSTOM",
  );

  return (
    <Menu isLazy>
      <MenuButton
        as={Button}
        size="sm"
        colorScheme="green"
        leftIcon={<FaPlus />}
        rightIcon={<FaChevronDown />}
        isLoading={isLoading}
      >
        Add step
      </MenuButton>
      <MenuList maxH="60vh" overflowY="auto" zIndex="dropdown">
        <MenuGroup title="Common">
          {common.map((type) => (
            <MenuItem key={type} onClick={() => onAdd(type)}>
              <StepTypeLabel stepType={type} />
            </MenuItem>
          ))}
        </MenuGroup>
        <MenuDivider />
        <MenuGroup title="More">
          {others.map((type) => (
            <MenuItem key={type} onClick={() => onAdd(type)}>
              <StepTypeLabel stepType={type} />
            </MenuItem>
          ))}
          <MenuItem onClick={() => onAdd("CUSTOM")}>
            <StepTypeLabel stepType="CUSTOM" />
          </MenuItem>
        </MenuGroup>
      </MenuList>
    </Menu>
  );
};

/**
 * The workup or purification card on the reaction page: the reaction's steps,
 * an "Add step" menu, and saving or using workbook templates.
 * @param {{
 *   workgroup: string,
 *   workbook: string,
 *   reaction: string,
 *   category: "WORKUP" | "PURIFICATION",
 * }} props
 */
export const WorkupCard = ({ workgroup, workbook, reaction, category }) => {
  const reactionNames = useMemo(
    () => ({ workgroup, workbook, reaction }),
    [workgroup, workbook, reaction],
  );
  const workbookNames = useMemo(
    () => ({ workgroup, workbook }),
    [workgroup, workbook],
  );
  const typeLabel = CATEGORY_LABELS[category];

  const { workups: api } = useBackendApi();
  const steps = useWorkupSteps(reactionNames, category);
  const definitions = useStepDefinitions();
  const templates = useWorkupTemplates(workbookNames, category);
  const templatesModal = useDisclosure();
  const saveModal = useDisclosure();
  const [isAdding, setIsAdding] = useState(false);
  const [error, setError] = useState(null);
  // success messages show in the card rather than as Chakra toasts: the reaction
  // page has several React roots, and each one would draw every toast
  const [notice, setNotice] = useState(null);

  useEffect(() => {
    if (!notice) return undefined;
    const timeout = setTimeout(() => setNotice(null), NOTICE_DURATION_MS);
    return () => clearTimeout(timeout);
  }, [notice]);

  if (steps.error || definitions.error) {
    return (
      <Alert status="error" fontSize="sm">
        <AlertIcon />
        Couldn't load the {typeLabel}. Try refreshing the page.
      </Alert>
    );
  }
  if (!steps.data || !definitions.data) {
    return <Spinner size="sm" aria-label={`Loading ${typeLabel}`} />;
  }

  const { canEdit } = steps.data;
  const stepList = steps.data.steps;

  /** Runs a write, showing its error in the card if it fails. */
  const run = async (action) => {
    setError(null);
    try {
      await action();
      return true;
    } catch (e) {
      console.error(e);
      setError(await describeWorkupError(e));
      return false;
    }
  };

  const handleAdd = (stepType) =>
    run(async () => {
      setIsAdding(true);
      try {
        await api.addStep(reactionNames, category, stepType);
        await steps.mutate();
      } finally {
        setIsAdding(false);
      }
    });

  const handleMove = (index, direction) =>
    run(async () => {
      const ids = stepList.map((step) => step.id);
      const target = index + direction;
      [ids[index], ids[target]] = [ids[target], ids[index]];
      const updated = await api.reorderSteps(reactionNames, category, ids);
      await steps.mutate(updated, { revalidate: false });
    });

  const handleDelete = (step) =>
    run(async () => {
      await api.deleteStep(reactionNames, step.id);
      await steps.mutate();
    });

  // keeps the cached list in step with a saved step, without refetching
  const handleSaved = (saved) =>
    steps.mutate(
      (current) =>
        current && {
          ...current,
          steps: current.steps.map((step) =>
            step.id === saved.id ? saved : step,
          ),
        },
      { revalidate: false },
    );

  const handleUseTemplate = async (template) => {
    const updated = await api.applyTemplate(reactionNames, template.id);
    await steps.mutate(updated, { revalidate: false });
    setNotice(`Added the steps from "${template.name}"`);
  };

  const handleDeleteTemplate = async (template) => {
    await api.deleteTemplate(workbookNames, template.id);
    await templates.mutate();
  };

  const handleSaveTemplate = async (values) => {
    await api.saveAsTemplate(reactionNames, category, values);
    await templates.mutate();
    setNotice("Template saved");
  };

  return (
    <Stack spacing={3}>
      <HStack spacing={2} flexWrap="wrap">
        {canEdit && (
          <AddStepMenu
            category={category}
            stepTypes={Object.keys(definitions.data)}
            onAdd={handleAdd}
            isLoading={isAdding}
          />
        )}
        <Button
          size="sm"
          variant="outline"
          leftIcon={<FaCopy />}
          onClick={templatesModal.onOpen}
        >
          {canEdit ? "Start from template" : "Templates"}
        </Button>
        {stepList.length > 0 && (
          <Button
            size="sm"
            variant="outline"
            leftIcon={<FaSave />}
            onClick={saveModal.onOpen}
          >
            Save as template
          </Button>
        )}
      </HStack>

      {notice && (
        <Alert status="success" fontSize="sm" role="status">
          <AlertIcon />
          {notice}
        </Alert>
      )}

      {error && (
        <Alert status="error" fontSize="sm">
          <AlertIcon />
          {error}
        </Alert>
      )}

      {stepList.length === 0 ? (
        <Box borderWidth="1px" borderStyle="dashed" borderRadius="md" p={4}>
          <Text fontSize="sm" color="gray.600">
            {canEdit
              ? `No ${typeLabel} steps yet. Use "Add step" or start from a template.`
              : `No ${typeLabel} steps recorded.`}
          </Text>
        </Box>
      ) : (
        <Stack spacing={2}>
          {stepList.map((step, index) => (
            <StepCard
              key={step.id}
              step={step}
              position={index + 1}
              isFirst={index === 0}
              isLast={index === stepList.length - 1}
              definition={definitions.data[step.stepType] ?? EMPTY_DEFINITION}
              canEdit={canEdit}
              reactionNames={reactionNames}
              onMove={(direction) => handleMove(index, direction)}
              onDelete={() => handleDelete(step)}
              onSaved={handleSaved}
            />
          ))}
        </Stack>
      )}

      <TemplatesModal
        isOpen={templatesModal.isOpen}
        onClose={templatesModal.onClose}
        templates={templates.data}
        isLoading={!templates.data && !templates.error}
        loadError={Boolean(templates.error)}
        canUse={canEdit}
        onUse={handleUseTemplate}
        onDelete={handleDeleteTemplate}
        typeLabel={typeLabel}
      />
      <SaveTemplateModal
        isOpen={saveModal.isOpen}
        onClose={saveModal.onClose}
        onSave={handleSaveTemplate}
        typeLabel={typeLabel}
      />
    </Stack>
  );
};
