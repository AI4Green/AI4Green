import {
  Alert,
  AlertIcon,
  Button,
  FormControl,
  FormErrorMessage,
  FormLabel,
  HStack,
  Input,
  Modal,
  ModalBody,
  ModalCloseButton,
  ModalContent,
  ModalFooter,
  ModalHeader,
  ModalOverlay,
  Spinner,
  Stack,
  Text,
  Textarea,
  Wrap,
  WrapItem,
} from "@chakra-ui/react";
import { Field, Form, Formik } from "formik";
import { useState } from "react";
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

/**
 * Asks for a name and description, then saves the card's steps as a template.
 * @param {{ isOpen: boolean, onClose: () => void, onSave: (values: object) => Promise<void>, typeLabel: string }} props
 */
export const SaveTemplateModal = ({ isOpen, onClose, onSave, typeLabel }) => {
  const [error, setError] = useState(null);

  const handleSubmit = async (values) => {
    setError(null);
    try {
      await onSave({
        name: values.name.trim(),
        description: values.description.trim(),
      });
      onClose();
    } catch (e) {
      console.error(e);
      setError(await describeWorkupError(e));
    }
  };

  return (
    <Modal isOpen={isOpen} onClose={onClose}>
      <ModalOverlay />
      <ModalContent>
        <Formik
          initialValues={{ name: "", description: "" }}
          validationSchema={saveTemplateSchema}
          onSubmit={handleSubmit}
        >
          {({ isSubmitting }) => (
            <Form noValidate>
              <ModalHeader>Save {typeLabel} as template</ModalHeader>
              <ModalCloseButton />
              <ModalBody>
                <Stack spacing={3}>
                  <Text fontSize="sm" color="gray.600">
                    Saves the steps in order, without the values, so anyone in
                    this workbook can start from them.
                  </Text>
                  <Field name="name">
                    {({ field, meta }) => (
                      <FormControl
                        isRequired
                        isInvalid={Boolean(meta.touched && meta.error)}
                      >
                        <FormLabel>Template name</FormLabel>
                        <Input {...field} maxLength={200} />
                        <FormErrorMessage>{meta.error}</FormErrorMessage>
                      </FormControl>
                    )}
                  </Field>
                  <Field name="description">
                    {({ field, meta }) => (
                      <FormControl
                        isInvalid={Boolean(meta.touched && meta.error)}
                      >
                        <FormLabel>Description</FormLabel>
                        <Textarea {...field} rows={3} maxLength={1000} />
                        <FormErrorMessage>{meta.error}</FormErrorMessage>
                      </FormControl>
                    )}
                  </Field>
                  {error && (
                    <Alert status="error" fontSize="sm">
                      <AlertIcon />
                      {error}
                    </Alert>
                  )}
                </Stack>
              </ModalBody>
              <ModalFooter>
                <Button variant="ghost" mr={2} onClick={onClose}>
                  Cancel
                </Button>
                <Button
                  colorScheme="green"
                  type="submit"
                  isLoading={isSubmitting}
                >
                  Save template
                </Button>
              </ModalFooter>
            </Form>
          )}
        </Formik>
      </ModalContent>
    </Modal>
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
      borderWidth="1px"
      borderRadius="md"
      p={3}
      spacing={2}
      aria-label={template.name}
    >
      <HStack align="start">
        <Stack spacing={0} flex="1">
          <Text fontWeight="semibold">{template.name}</Text>
          {template.description && (
            <Text fontSize="sm" color="gray.600">
              {template.description}
            </Text>
          )}
          {template.creator && (
            <Text fontSize="xs" color="gray.500">
              Saved by {template.creator}
            </Text>
          )}
        </Stack>
        {canUse && (
          <Button
            size="sm"
            colorScheme="green"
            onClick={() => run(onUse, "use")}
            isLoading={busy === "use"}
            isDisabled={busy !== null}
          >
            Use
          </Button>
        )}
        {template.canDelete && (
          <Button
            size="sm"
            variant="ghost"
            colorScheme="red"
            onClick={() => run(onDelete, "delete")}
            isLoading={busy === "delete"}
            isDisabled={busy !== null}
          >
            Delete
          </Button>
        )}
      </HStack>
      <Wrap spacing={2} fontSize="xs" color="gray.700">
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
    onClose();
  });

  let body;
  if (loadError) {
    body = <Text>Couldn't load the templates. Try refreshing the page.</Text>;
  } else if (isLoading) {
    body = <Spinner />;
  } else if (!templates?.length) {
    body = (
      <Text color="gray.600">
        No {typeLabel} templates in this workbook yet. Add some steps, then use
        "Save as template".
      </Text>
    );
  } else {
    body = (
      <Stack as="ul" spacing={3} m={0} p={0}>
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
    <Modal isOpen={isOpen} onClose={onClose} size="lg" scrollBehavior="inside">
      <ModalOverlay />
      <ModalContent>
        <ModalHeader textTransform="capitalize">
          {typeLabel} templates
        </ModalHeader>
        <ModalCloseButton />
        <ModalBody>
          <Stack spacing={3}>
            {canUse && templates?.length > 0 && (
              <Text fontSize="sm" color="gray.600">
                Using a template adds its steps after any you already have.
              </Text>
            )}
            {error && (
              <Alert status="error" fontSize="sm">
                <AlertIcon />
                {error}
              </Alert>
            )}
            {body}
          </Stack>
        </ModalBody>
        <ModalFooter>
          <Button onClick={onClose}>Close</Button>
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
};
