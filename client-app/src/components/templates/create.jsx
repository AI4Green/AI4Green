import {
  Alert,
  AlertIcon,
  Box,
  Button,
  CloseButton,
  Heading,
  HStack,
  useToast,
} from "@chakra-ui/react";
import { TemplateSelectorForm } from "components/templates";

export const CreateFromTemplateInline = ({
  title,
  templateLabel,
  formRef,
  templates,
  initialValues,
  feedback,
  isLoading,
  onSubmit,
  onClose,
  onCreateTemplate,
  createLabel = "Create",
  isOpen,
}) => {
  if (!isOpen) return null;

  return (
    <Box borderWidth="1px" borderRadius="md" p={4} bg="white">
      <HStack justify="space-between" mb={3}>
        <Heading size="sm">{title}</Heading>

        {onClose && <CloseButton onClick={onClose} />}
      </HStack>

      <TemplateSelectorForm
        formRef={formRef}
        templateLabel={templateLabel}
        templates={templates}
        initialValues={initialValues}
        feedback={feedback}
        onSubmit={onSubmit}
      />

      <HStack mt={4}>
        <Button
          colorScheme="green"
          onClick={() => formRef.current?.handleSubmit()}
          isLoading={isLoading}
        >
          {createLabel}
        </Button>

        {onCreateTemplate && (
          <Button variant="outline" onClick={onCreateTemplate}>
            Create new template
          </Button>
        )}
      </HStack>
    </Box>
  );
};
