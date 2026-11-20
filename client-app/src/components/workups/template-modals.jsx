import { HStack, Spinner, Stack, Text, Wrap, WrapItem } from "@chakra-ui/react";
import { TextAreaField, TextField } from "components/inventory/form-fields";
import {
  ErrorMessage,
  INV,
  InvButton,
  InvModal,
  Muted,
  underline,
} from "components/inventory/ui";
import { Form, Formik } from "formik";
import { useId, useState } from "react";
import * as Yup from "yup";

import { describeWorkupError } from "./errors";
import { StepTypeLabel } from "./step-type-label";

const saveTemplateSchema = Yup.object({
  name: Yup.string()
    .trim()
    .required("Give the template a name.")
    .max(200, "Use 200 characters or fewer."),
  description: Yup.string().max(1000, "Use 1000 characters or fewer."),
});

/** "workup" -> "Workup", for headings. */
const capitalise = (text) => text.charAt(0).toUpperCase() + text.slice(1);

/**
 * Asks for a name and description, then saves the card's steps as a template.
 * @param {{ isOpen: boolean, onClose: () => void, onSave: (values: object) => Promise<void>, typeLabel: string }} props
 */
export const SaveTemplateModal = ({ isOpen, onClose, onSave, typeLabel }) => {
  const formId = `save-template-${useId()}`;
  const [error, setError] = useState(null);
  const [isSaving, setIsSaving] = useState(false);

  const handleClose = () => {
    setError(null);
    onClose();
  };

  const handleSubmit = async (values) => {
    setError(null);
    setIsSaving(true);
    try {
      await onSave({
        name: values.name.trim(),
        description: values.description.trim(),
      });
      handleClose();
    } catch (e) {
      console.error(e);
      setError(await describeWorkupError(e));
    } finally {
      setIsSaving(false);
    }
  };

  return (
    <InvModal
      isOpen={isOpen}
      onClose={handleClose}
      title={`Save ${typeLabel} as template`}
      actionLabel="Save template"
      formId={formId}
      isLoading={isSaving}
    >
      {error && <ErrorMessage>{error}</ErrorMessage>}
      <Formik
        initialValues={{ name: "", description: "" }}
        validationSchema={saveTemplateSchema}
        onSubmit={handleSubmit}
      >
        <Stack as={Form} id={formId} noValidate spacing="16px">
          <Text m={0} color={INV.muted}>
            Saves the steps in order, without the values, so anyone in this
            workbook can start from them.
          </Text>
          <TextField
            name="name"
            label="Template name"
            maxLength={200}
            isRequired
            autoFocus
          />
          <TextAreaField
            name="description"
            label="Description"
            rows={3}
            maxLength={1000}
          />
        </Stack>
      </Formik>
    </InvModal>
  );
};

/** One template in the list, with its steps and buttons. */
const TemplateRow = ({ template, canUse, onUse, onDelete }) => {
  const [busy, setBusy] = useState(null);

  const run = async (action, name) => {
    setBusy(name);
    try {
      await action(template);
    } finally {
      setBusy(null);
    }
  };

  return (
    <Stack
      as="li"
      listStyleType="none"
      border="1px solid"
      borderColor={INV.line}
      borderRadius="12px"
      p="12px"
      spacing="8px"
      aria-label={template.name}
    >
      <HStack align="start" spacing="8px">
        <Stack spacing="2px" flex="1">
          <Text m={0} fontWeight={600}>
            {template.name}
          </Text>
          {template.description && (
            <Text m={0} color={INV.muted}>
              {template.description}
            </Text>
          )}
          {template.creator && <Muted>Saved by {template.creator}</Muted>}
        </Stack>
        {canUse && (
          <InvButton
            primary
            onClick={() => run(onUse, "use")}
            isLoading={busy === "use"}
            isDisabled={busy !== null}
          >
            Use
          </InvButton>
        )}
        {template.canDelete && (
          <InvButton
            color={INV.danger}
            borderColor={INV.danger}
            _hover={{ bg: INV.dangerSoft, ...underline }}
            onClick={() => run(onDelete, "delete")}
            isLoading={busy === "delete"}
            isDisabled={busy !== null}
          >
            Delete
          </InvButton>
        )}
      </HStack>
      <Wrap spacing="8px" fontSize="12px" color={INV.muted}>
        {template.steps.map((stepType, index) => (
          <WrapItem key={`${stepType}-${index}`}>
            <StepTypeLabel stepType={stepType} />
          </WrapItem>
        ))}
      </Wrap>
    </Stack>
  );
};

/**
 * Lists the workbook's templates for this card. "Use" appends the template's
 * steps to the reaction; a template's creator can also delete it.
 * @param {{
 *   isOpen: boolean,
 *   onClose: () => void,
 *   templates: object[] | undefined,
 *   isLoading: boolean,
 *   loadError: boolean,
 *   canUse: boolean,
 *   onUse: (template: object) => Promise<void>,
 *   onDelete: (template: object) => Promise<void>,
 *   typeLabel: string,
 * }} props
 */
export const TemplatesModal = ({
  isOpen,
  onClose,
  templates,
  isLoading,
  loadError,
  canUse,
  onUse,
  onDelete,
  typeLabel,
}) => {
  const [error, setError] = useState(null);

  const handleClose = () => {
    setError(null);
    onClose();
  };

  const wrap = (action) => async (template) => {
    setError(null);
    try {
      await action(template);
    } catch (e) {
      console.error(e);
      setError(await describeWorkupError(e));
    }
  };

  const handleUse = wrap(async (template) => {
    await onUse(template);
    handleClose();
  });

  let body;
  if (loadError) {
    body = (
      <ErrorMessage mb={0}>
        Couldn't load the templates. Try refreshing the page.
      </ErrorMessage>
    );
  } else if (isLoading) {
    body = <Spinner color={INV.accent} aria-label="Loading templates" />;
  } else if (!templates?.length) {
    body = (
      <Text m={0} color={INV.muted}>
        No {typeLabel} templates in this workbook yet. Add some steps, then use
        "Save as template".
      </Text>
    );
  } else {
    body = (
      <Stack as="ul" spacing="12px" m={0} p={0}>
        {templates.map((template) => (
          <TemplateRow
            key={template.id}
            template={template}
            canUse={canUse}
            onUse={handleUse}
            onDelete={wrap(onDelete)}
          />
        ))}
      </Stack>
    );
  }

  return (
    <InvModal
      isOpen={isOpen}
      onClose={handleClose}
      title={`${capitalise(typeLabel)} templates`}
      size="lg"
    >
      <Stack spacing="12px">
        {canUse && templates?.length > 0 && (
          <Text m={0} color={INV.muted}>
            Using a template adds its steps after any you already have.
          </Text>
        )}
        {error && <ErrorMessage mb={0}>{error}</ErrorMessage>}
        {body}
      </Stack>
    </InvModal>
  );
};
