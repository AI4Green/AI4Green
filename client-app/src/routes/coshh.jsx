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
import { Formik, Form } from "formik";
import { useBackendApi } from "contexts";
import { useProjectTypesList } from "api";
import { FormikInput, MultiSelectField } from "components/core/forms";

export const CoshhCreateModal = () => {
  const { reactionId } = useParams();

  const navigate = useNavigate();
  const location = useLocation();

  const { projects: action } = useBackendApi();
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

      const response = await action.create({
        reaction_id: Number(reactionId),
        template_id: Number(values.templateId[0]),
      });

      if (response?.status === 200 || response?.status === 201) {
        toast({
          title: "COSHH form created",
          status: "success",
          duration: GLOBAL_PARAMETERS.ToastDuration,
          isClosable: true,
          position: "top",
        });

        handleReset();

        navigate(`/coshh/form/${response.data.id}/edit`);
      }
    } catch (e) {
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
  const { coshh: api } = useBackendApi();
  const [loading, setLoading] = useState(false);
  const [form, setForm] = useState(null);

  return (
    <VStack align="stretch" spacing={4} p={5}>
      <Text fontSize="xl" fontWeight="bold">
        Editing COSHH Form: {formId}
      </Text>
      {/* Your form fields would go here */}
      <Button colorScheme="green">Save Changes</Button>
    </VStack>
  );
};

export const COSHH = () => {
  return (
    <Routes>
      <Route path="new/:reactionId" element={<CoshhCreateModal />} />

      {/* When the URL is /coshh/form/:formId,
         the Modal is gone because this route doesn't render it.
      */}
      <Route path="form/:formId" element={<CoshhForm />} />
    </Routes>
  );
};
