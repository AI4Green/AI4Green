import React, { useRef, useState } from "react";
import {
  Box,
  Button,
  CloseButton,
  Heading,
  HStack,
  useToast,
} from "@chakra-ui/react";

import { useBackendApi } from "contexts";
import { useWorkupTemplatesList } from "api/templates";

import { WorkupCreateFormContent } from "./form";

export const WorkupCreateInline = ({
  reactionId,
  workgroupName,
  workbookName,
  onCreated,
  onCreateTemplate,
}) => {
  const {
    formRef,
    templates,
    initialValues,
    feedback,
    isLoading,
    handleSubmit,
  } = useWorkupCreateForm({
    reactionId,
    workgroupName,
    workbookName,
    onCreated,
  });

  return (
    <Box borderWidth="1px" borderRadius="md" p={4} bg="white">
      <Heading size="sm" mb={3}>
        Create Workup
      </Heading>

      <WorkupCreateFormContent
        formRef={formRef}
        templates={templates}
        initialValues={initialValues}
        feedback={feedback}
        onSubmit={handleSubmit}
      />

      <HStack mt={4}>
        <Button
          colorScheme="green"
          onClick={() => formRef.current?.handleSubmit()}
          isLoading={isLoading}
        >
          Create
        </Button>

        <Button variant="outline" onClick={onCreateTemplate}>
          Create new template
        </Button>
      </HStack>
    </Box>
  );
};
