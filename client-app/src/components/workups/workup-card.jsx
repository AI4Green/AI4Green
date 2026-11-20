import {
  Box,
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
import {
  EmptyState,
  ErrorMessage,
  INV,
  InvButton,
  underline,
} from "components/inventory/ui";
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
// menu entries get the same hover/focus underline as the buttons (WCAG 1.4.1)
const menuItemProps = {
  fontSize: "14px",
  color: INV.ink,
  _hover: { bg: INV.accentSoft, ...underline },
  _focus: { bg: INV.accentSoft, ...underline },
};

/** The "Add step" menu: this card's common step types first, then the rest. */
const AddStepMenu = ({ category, stepTypes, onAdd, isLoading }) => {
  const common = COMMON_STEP_TYPES[category] ?? [];
  const others = stepTypes.filter(
    (type) => !common.includes(type) && type !== "CUSTOM",
  );

  return (
    <Menu isLazy>
      <MenuButton
        as={InvButton}
        primary
        leftIcon={<FaPlus />}
        rightIcon={<FaChevronDown />}
        isLoading={isLoading}
      >
        Add step
      </MenuButton>
      <MenuList
        maxH="60vh"
        overflowY="auto"
        zIndex="dropdown"
        borderRadius="12px"
        borderColor={INV.line}
      >
        <MenuGroup title="Common" color={INV.muted} fontSize="12px">
          {common.map((type) => (
            <MenuItem key={type} onClick={() => onAdd(type)} {...menuItemProps}>
              <StepTypeLabel stepType={type} />
            </MenuItem>
          ))}
        </MenuGroup>
        <MenuDivider borderColor={INV.line} />
        <MenuGroup title="More" color={INV.muted} fontSize="12px">
          {others.map((type) => (
            <MenuItem key={type} onClick={() => onAdd(type)} {...menuItemProps}>
              <StepTypeLabel stepType={type} />
            </MenuItem>
          ))}
          <MenuItem onClick={() => onAdd("CUSTOM")} {...menuItemProps}>
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
  // the step the user just added, which opens ready to fill in; steps loaded
  // with the page stay collapsed
  const [addedStepId, setAddedStepId] = useState(null);
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
      <ErrorMessage mb={0}>
        Couldn't load the {typeLabel}. Try refreshing the page.
      </ErrorMessage>
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
        const added = await api.addStep(reactionNames, category, stepType);
        setAddedStepId(added.id);
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
        <InvButton leftIcon={<FaCopy />} onClick={templatesModal.onOpen}>
          {canEdit ? "Start from template" : "Templates"}
        </InvButton>
        {stepList.length > 0 && (
          <InvButton leftIcon={<FaSave />} onClick={saveModal.onOpen}>
            Save as template
          </InvButton>
        )}
      </HStack>

      {notice && (
        <Text
          role="status"
          m={0}
          px="12px"
          py="10px"
          borderRadius="8px"
          bg={INV.accentSoft}
          color={INV.accentHover}
          fontSize="14px"
        >
          {notice}
        </Text>
      )}

      {error && <ErrorMessage mb={0}>{error}</ErrorMessage>}

      {stepList.length === 0 ? (
        <Box
          border="1px dashed"
          borderColor={INV.lineStrong}
          borderRadius="12px"
          bg="white"
        >
          <EmptyState py="20px" fontSize="14px">
            {canEdit
              ? `No ${typeLabel} steps yet. Use "Add step" or start from a template.`
              : `No ${typeLabel} steps recorded.`}
          </EmptyState>
        </Box>
      ) : (
        <Stack spacing="8px">
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
              defaultExpanded={step.id === addedStepId}
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
