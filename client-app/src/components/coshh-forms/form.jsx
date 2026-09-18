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
import { useCreateFormFromTemplate } from "components/templates";
import { TemplateSelectorForm } from "components/templates/select.jsx";
import { TemplateCreateModalInline } from "components/templates/create.jsx";
import { CoshhCreateInline } from "components/coshh-forms";

// modal for routed/standalone page
export const CoshhCreateModal = ({
  reactionId,
  onCreated,
  isOpen,
  onClose,
}) => {
  const { formRef, initialValues, feedback, isLoading, handleSubmit } =
    useCreateFormFromTemplate({ templateType: "COSHH", reactionId, onCreated });

  const { data: templates = [] } = useCoshhFormsList();

  const body = (
    <TemplateSelectorForm
      formRef={formRef}
      label={"COSHH"}
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
    isOwner: true, // todo: control non creator functionality
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
