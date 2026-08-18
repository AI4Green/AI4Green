import {
  Alert,
  AlertIcon,
  HStack,
  Icon,
  Text,
  useToast,
  VStack,
} from "@chakra-ui/react";
import { useProjectTypesList } from "api";
import {
  FormikInput,
  MultiSelectField,
  TextAreaField,
} from "components/core/forms";
import { Modal, useModalState } from "components/core/modal";
import { GLOBAL_PARAMETERS, STAGES } from "constants";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { useEffect, useRef } from "react";
import { useTranslation } from "react-i18next";
import { FaLayerGroup } from "react-icons/fa";
import { useLocation, useNavigate, useSearchParams } from "react-router-dom";
import { object, string } from "yup";

export const CreateOrEditProjectTypeModal = () => {
  const [searchParams] = useSearchParams();
  const id = searchParams.get("id");
  const isEditAction = searchParams.get("action") === "edit";

  const navigate = useNavigate();
  const location = useLocation();

  const { projectTypes: action } = useBackendApi();
  const { data: projectTypes, mutate } = useProjectTypesList();

  const { t } = useTranslation();
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

  const projectType = isEditAction
    ? projectTypes?.find((projectType) => projectType.id === Number(id))
    : null;

  useEffect(() => {
    setIsModalOpen(true);
  }, [id, isEditAction, setIsModalOpen]);

  const initialValues = () => {
    return projectType
      ? {
          name: projectType.name,
          description: projectType.description,
          source: [],
        }
      : {
          name: "",
          description: "",
          source: [],
        };
  };

  const handleSubmit = async (values) => {
    try {
      setIsLoading(true);
      const model = {
        name: values.name,
        description: values.description,
        sourceId: values.source.length > 0 ? Number(values.source[0]) : null,
      };
      const response = !projectType
        ? await action.create({ values: model })
        : await action.edit({ values: model, id: projectType.id });
      setIsLoading(false);

      if (response && (response.status === 204 || response.status === 200)) {
        toast({
          title: `Project Type ${!projectType ? "created" : "updated"}`,
          status: "success",
          duration: GLOBAL_PARAMETERS.ToastDuration,
          isClosable: true,
          position: "top",
        });
        handleReset();
        await mutate();
      }
    } catch (e) {
      switch (e?.response?.status) {
        case 400: {
          setFeedback({
            status: "error",
            message: t("feedback.error_400"),
          });
          break;
        }
        default: {
          setFeedback({
            status: "error",
            message: t("feedback.error"),
          });
          break;
        }
      }
    }
  };

  const modalBody = (
    <Formik
      enableReinitialize
      innerRef={formRef}
      initialValues={initialValues()}
      onSubmit={handleSubmit}
      validationSchema={validationSchema(projectTypes, projectType)}
    >
      <Form noValidate>
        <VStack align="stretch" spacing={4}>
          {feedback && (
            <Alert status={feedback.status}>
              <AlertIcon />
              {feedback.message}
            </Alert>
          )}
          <HStack spacing={5} align="start">
            <Icon
              as={FaLayerGroup}
              color={projectType ? "blue.500" : "green.500"}
              fontSize="5xl"
            />
            <VStack w="full" align="start">
              <FormikInput name="name" label="Project type name" isRequired />
              <TextAreaField
                name="description"
                label="Project type description"
              />
              <Text fontSize="xs" color="gray.500" as="i">
                Optionally import sections and fields from an existing project
                type. Leave blank to create an empty project type.
              </Text>
              <MultiSelectField
                name="source"
                label="Import from"
                options={projectTypes
                  .filter((x) => x.stage === STAGES.Ready)
                  .map((projectType) => ({
                    label: projectType.name,
                    value: String(projectType.id),
                    description: projectType.description,
                  }))}
              />
            </VStack>
          </HStack>
        </VStack>
      </Form>
    </Formik>
  );

  return (
    <Modal
      body={modalBody}
      title={`${!isEditAction ? "Create" : "Edit"} Project Type`}
      actionBtnCaption={!isEditAction ? "Create" : "Update"}
      onAction={() => formRef.current.handleSubmit()}
      actionBtnColorScheme={!isEditAction ? "green" : "blue"}
      isLoading={isLoading}
      isOpen={isModalOpen}
      onClose={handleReset}
    />
  );
};

const validationSchema = (projectTypes, current) => {
  const existingNames = projectTypes
    .filter((x) => !current || x.id !== current.id)
    .map((projectType) => projectType.name);

  return object().shape({
    name: string()
      .notOneOf(existingNames, "Project type name already exist")
      .required("Project name required"),
    description: string().required("Project type description required"),
  });
};
