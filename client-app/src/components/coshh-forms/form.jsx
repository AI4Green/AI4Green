import { useNavigate, useParams, useLocation } from "react-router-dom";
import { Modal, useModalState } from "components/core/modal";
import { useRef, useEffect, useState } from "react";
import {
  VStack,
  HStack,
  Box,
  Heading,
  Button,
  CloseButton,
  Alert,
  AlertIcon,
  useToast,
} from "@chakra-ui/react";
import { SectionForm } from "components/section-form";
import { Formik, Form } from "formik";
import { useBackendApi } from "contexts";
import { useCoshhFormsList, useProject } from "api";
import { FormikInput, MultiSelectField } from "components/core/forms";

const useCoshhCreateForm = ({
  reactionId,
  workgroupName,
  workbookName,
  onCreated,
}) => {
  const { projects: api } = useBackendApi();
  const { data: templates = [] } = useCoshhFormsList();
  const toast = useToast();
  const formRef = useRef();

  const [isLoading, setIsLoading] = useState(false);
  const [feedback, setFeedback] = useState(null);

  const initialValues = { templateId: [] };

  const handleSubmit = async (values) => {
    try {
      setIsLoading(true);
      setFeedback(null);

      const response = await api.create({
        reactionId,
        workgroupName,
        workbookName,
        templateId: Number(values.templateId[0]),
        templateType: "COSHH",
      });

      if (!response.ok) throw new Error("Failed to create COSHH form");

      const data = await response.json();

      if (data?.id || data?.uuid) {
        toast({
          title: "COSHH form created",
          status: "success",
          duration: 10000,
          isClosable: true,
          position: "top",
        });
        onCreated(data);
      }
    } catch (e) {
      console.error(e);
      setFeedback({ status: "error", message: "Failed to create COSHH form" });
    } finally {
      setIsLoading(false);
    }
  };

  return {
    formRef,
    templates,
    initialValues,
    feedback,
    isLoading,
    handleSubmit,
  };
};

const CoshhCreateFormContent = ({
  formRef,
  templates,
  initialValues,
  feedback,
  onSubmit,
}) => (
  <Formik innerRef={formRef} initialValues={initialValues} onSubmit={onSubmit}>
    {() => (
      <Form noValidate>
        <VStack spacing={4} align="stretch">
          {feedback && (
            <Alert status={feedback.status}>
              <AlertIcon />
              {feedback.message}
            </Alert>
          )}
          <MultiSelectField
            isRequired
            name="templateId"
            label="COSHH Template"
            options={templates.map((template) => ({
              label: template.name,
              value: String(template.id),
              description: template.description,
            }))}
          />
        </VStack>
      </Form>
    )}
  </Formik>
);

// modal for routed/standalone page
export const CoshhCreateModal = ({
  reactionId,
  onCreated,
  isOpen,
  onClose,
}) => {
  const {
    formRef,
    templates,
    initialValues,
    feedback,
    isLoading,
    handleSubmit,
  } = useCoshhCreateForm({ reactionId, onCreated });

  const body = (
    <CoshhCreateFormContent
      formRef={formRef}
      templates={templates}
      initialValues={initialValues}
      feedback={feedback}
      onSubmit={handleSubmit}
    />
  );

  return (
    <Modal
      body={body}
      title="Create COSHH Form"
      actionBtnCaption="Create"
      actionBtnColorScheme="green"
      onAction={() => formRef.current?.handleSubmit()}
      isLoading={isLoading}
      isOpen={isOpen}
      onClose={onClose}
    />
  );
};

// for embedding into reaction constructor
export const CoshhCreateInline = ({
  reactionId,
  workgroupName,
  workbookName,
  onCreated,
  isOpen,
  onClose,
}) => {
  const {
    formRef,
    templates,
    initialValues,
    feedback,
    isLoading,
    handleSubmit,
  } = useCoshhCreateForm({
    reactionId,
    workgroupName,
    workbookName,
    onCreated,
  });

  if (!isOpen) return null;

  return (
    <Box borderWidth="1px" borderRadius="md" p={4} bg="white">
      <HStack justify="space-between" mb={3}>
        <Heading size="sm">Create COSHH Form</Heading>
        <CloseButton onClick={onClose} />
      </HStack>

      <CoshhCreateFormContent
        formRef={formRef}
        templates={templates}
        initialValues={initialValues}
        feedback={feedback}
        onSubmit={handleSubmit}
      />

      <Button
        mt={4}
        colorScheme="green"
        onClick={() => formRef.current?.handleSubmit()}
        isLoading={isLoading}
      >
        Create
      </Button>
    </Box>
  );
};

// routed standalone page
export const RoutedCoshhCreateModal = () => {
  const { reactionId } = useParams();
  const navigate = useNavigate();
  const location = useLocation();
  const formRef = useRef();

  const { isModalOpen, setIsModalOpen, handleReset } = useModalState(
    location,
    navigate,
    formRef,
  );

  useEffect(() => {
    setIsModalOpen(true);
  }, [setIsModalOpen]);

  return (
    <CoshhCreateModal
      reactionId={reactionId}
      isOpen={isModalOpen}
      onClose={handleReset}
      onCreated={(data) => navigate(`/coshh/form/${data.id}/edit`)}
    />
  );
};

// embedded in the card
export const EmbeddedCoshh = ({
  reactionId,
  initialFormId = null,
  workbookName,
  workgroupName,
}) => {
  const [formId, setFormId] = useState(initialFormId);
  const [isCreateOpen, setIsCreateOpen] = useState(!initialFormId);

  if (formId) return <CoshhForm formId={formId} />;

  return (
    <CoshhCreateInline
      reactionId={reactionId}
      workgroupName={workgroupName}
      workbookName={workbookName}
      isOpen={isCreateOpen}
      onClose={() => setIsCreateOpen(false)}
      onCreated={(data) => {
        setFormId(data.id);
        setIsCreateOpen(false);
      }}
    />
  );
};

export const CoshhForm = ({ formId }) => {
  const { coshhForms: coshhFormsApi, projects: projectApi } = useBackendApi();
  const toast = useToast();

  const [data, setData] = useState(null);
  const [loading, setLoading] = useState(true);

  const { data: template } = useProject(formId);

  useEffect(() => {
    const loadCoshhData = async () => {
      try {
        setLoading(true);
        setData(template);
      } catch (e) {
        console.log(e);
        toast({ title: "Error loading form", status: "error" });
      } finally {
        setLoading(false);
      }
    };
    loadCoshhData();
  }, [formId, coshhFormsApi, toast]);

  const itemContext = {
    id: template.id,
    isOwner: true, // Usually true if they just created it
    type: "COSHH",
    approvalStatus: {
      stage: template.approvalStatus,
      permissions: "OwnerCanEdit",
    } || {
      // handle permissions better
      permissions: ["OWNER_CAN_EDIT"],
    },
    action: {
      save: async (formData) => projectApi.putForm(template.id, formData),
      mutate: () => window.location.reload(), // Or a more elegant SWR mutate
    },
  };
  return (
    <>
      {template.sections
        ?.sort((a, b) => a.sortOrder - b.sortOrder)
        .map((section) => (
          <SectionForm
            key={section.id}
            item={itemContext}
            form={section}
            isInstructor={false}
            breadcrumbItems={[
              {
                label: template.workgroup,
                href: `/workgroup/${template.workgroup}`,
                external: true,
              },
              {
                label: template.workbook,
                href: `/workgroup/${template.workgroup}`,
                external: true,
              },
              {
                label: template.reactionCode,
                href: `/sketcher/${template.workgroup}/${template.workbook}/${template.reactionCode}/no`,
                external: true,
              },
              { label: "COSHH", active: true },
            ]}
            headerItems={{
              title: section.name || "COSHH Assessment",
              subtitle: `Editing instance ${template.uuid}`,
              name: section.title || "SECTION",
            }}
          />
        ))}
    </>
  );
};
