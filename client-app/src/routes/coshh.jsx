import {
  Routes,
  Route,
  useNavigate,
  useParams,
  useLocation,
} from "react-router-dom";
import { Modal, useModalState } from "components/core/modal";
import { useRef, useEffect, useState } from "react";
import {
  VStack,
  Button,
  Text,
  Alert,
  AlertIcon,
  useToast,
} from "@chakra-ui/react";
import { SectionForm } from "components/section-form";
import { Formik, Form } from "formik";
import { useBackendApi } from "contexts";
import { useProjectTypesList, useProject } from "api";
import { FormikInput, MultiSelectField } from "components/core/forms";

export const CoshhCreateModal = () => {
  const { reactionId } = useParams();

  const navigate = useNavigate();
  const location = useLocation();

  const { projects: api } = useBackendApi();
  const { data: templates = [] } = useProjectTypesList();

  const toast = useToast();
  const formRef = useRef();

  const {
    isModalOpen,
    setIsModalOpen,
    isLoading,
    setIsLoading,
    feedback,
    setFeedback,
    handleReset,
  } = useModalState(location, navigate, formRef);

  useEffect(() => {
    setIsModalOpen(true);
  }, [setIsModalOpen]);

  const initialValues = {
    templateId: [],
  };

  const handleSubmit = async (values) => {
    try {
      setIsLoading(true);

      const response = await api.create({
        reactionId: reactionId,
        templateId: Number(values.templateId[0]),
        templateType: "COSHH", // best to send this here or have a dedicated route?
      });
      const data = await response.json();

      if (data?.id || data?.uuid) {
        toast({
          title: "COSHH form created",
          status: "success",
          duration: 10,
          isClosable: true,
          position: "top",
        });

        // handleReset();

        navigate(`/coshh/form/${data.id}/edit`);
      }
    } catch (e) {
      console.log(e);
      setFeedback({
        status: "error",
        message: "Failed to create COSHH form",
      });
    } finally {
      setIsLoading(false);
    }
  };

  const modalBody = (
    <Formik
      innerRef={formRef}
      initialValues={initialValues}
      onSubmit={handleSubmit}
    >
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

  return (
    <Modal
      body={modalBody}
      title="Create COSHH Form"
      actionBtnCaption="Create"
      actionBtnColorScheme="green"
      onAction={() => formRef.current.handleSubmit()}
      isLoading={isLoading}
      isOpen={isModalOpen}
      onClose={handleReset}
    />
  );
};

export const CoshhForm = () => {
  const { formId } = useParams();
  const { projectType: projectTypesApi, projects: projectApi } =
    useBackendApi();
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
  }, [formId, projectTypesApi, toast]);

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
              },
              {
                label: template.workbook,
                href: `/workgroup/${template.workgroup}`,
              },
              {
                label: template.reactionCode,
                href: `/sketcher/${template.workgroup}/${template.workbook}/${template.reationCode}/no`,
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

export const COSHH = () => {
  return (
    <Routes>
      <Route path="new/:reactionId" element={<CoshhCreateModal />} />

      {/* When the URL is /coshh/form/:formId,
         the Modal is gone because this route doesn't render it.
      */}
      <Route path="form/:formId/edit" element={<CoshhForm />} />
    </Routes>
  );
};
