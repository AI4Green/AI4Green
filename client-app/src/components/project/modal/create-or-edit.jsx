import {
  Alert,
  AlertIcon,
  HStack,
  Icon,
  useToast,
  VStack,
} from "@chakra-ui/react";
import { useProjectTypesList } from "api";
import { Badge } from "components/core/Badge";
import { FormikInput, MultiSelectField } from "components/core/forms";
import { Modal, useModalState } from "components/core/modal";
import { GLOBAL_PARAMETERS, STAGES, TITLE_ICON_COMPONENTS } from "constants";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { useEffect, useRef } from "react";
import { useTranslation } from "react-i18next";
import { useLocation, useNavigate, useSearchParams } from "react-router-dom";

import { validationSchema } from "./validation";

export const CreateOrEditModal = ({ list, mutate }) => {
  const [searchParams] = useSearchParams();
  const id = searchParams.get("id");
  const isEditAction = searchParams.get("action") === "edit";

  const navigate = useNavigate();
  const location = useLocation();

  const { projects: action } = useBackendApi();
  const { data: projectTypes } = useProjectTypesList();

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

  const project = isEditAction
    ? list?.find((project) => project.id === Number(id))
    : null;

  useEffect(() => {
    setIsModalOpen(true);
  }, [id, isEditAction, setIsModalOpen]);

  const initialValues = () => {
    return project
      ? {
          name: project.name,
          projectTypeId: [String(project.projectType.id)],
        }
      : {
          name: "",
          projectTypeId: [],
        };
  };

  const handleSubmit = async (values) => {
    try {
      setIsLoading(true);
      const model = {
        name: values.name,
        projectTypeId: Number(values.projectTypeId[0]),
      };
      const response = !project
        ? await action.create(model)
        : await action.edit(project.id, model);
      setIsLoading(false);

      if (response && (response.status === 204 || response.status === 200)) {
        toast({
          title: `Project ${!project ? "created" : "updated"}`,
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
      validationSchema={validationSchema(list, projectTypes)}
    >
      {({ values, setFieldValue }) => {
        return (
          <Form noValidate>
            <DeadlinesManager values={values} setFieldValue={setFieldValue} />
            <HStack spacing={4}>
              <VStack>
                <Icon
                  as={TITLE_ICON_COMPONENTS.Project}
                  color={isEditAction ? "blue.500" : "green.500"}
                  fontSize="5xl"
                />
                <Badge
                  colorScheme={isEditAction ? "blue" : "green"}
                  label="Project"
                  variant="outline"
                  fontSize="xxs"
                />
              </VStack>
              <VStack align="flex-start" flex={1}>
                {feedback && (
                  <Alert status={feedback.status}>
                    <AlertIcon />
                    {feedback.message}
                  </Alert>
                )}

                <VStack w="full" spacing={4}>
                  <FormikInput name="name" label="Project name" isRequired />
                  <MultiSelectField
                    isRequired
                    name="projectTypeId"
                    label="Project type"
                    options={projectTypes
                      .filter(
                        (projectType) => projectType.stage === STAGES.Ready,
                      )
                      .map((projectType) => ({
                        label: projectType.name,
                        value: String(projectType.id),
                        description: projectType.description,
                      }))}
                    isDisabled={!!project}
                  />
                </VStack>
              </VStack>
            </HStack>
          </Form>
        );
      }}
    </Formik>
  );

  return (
    <Modal
      body={modalBody}
      title={`${!project ? "Create" : "Edit"} Project`}
      actionBtnCaption={!project ? "Create" : "Update"}
      onAction={() => formRef.current.handleSubmit()}
      actionBtnColorScheme={!project ? "green" : "blue"}
      isLoading={isLoading}
      isOpen={isModalOpen}
      onClose={handleReset}
    />
  );
};

const DeadlinesManager = ({ values, setFieldValue }) => {
  useDeadline("planningDeadline", "startDate", 14, values, setFieldValue);
  useDeadline(
    "experimentDeadline",
    "planningDeadline",
    28,
    values,
    setFieldValue,
  );
  return null;
};

/**
 * Hook to set the deadline field value based on the baseField value
 * @param {*} field - field to be set
 * @param {*} baseField - field to be used as base for calculation
 * @param {*} daysToAdd - number of days to add to the baseField
 * @param {*} values - formik values
 * @param {*} setFieldValue - formik setFieldValue
 */
const useDeadline = (field, baseField, daysToAdd, values, setFieldValue) => {
  useEffect(() => {
    const deadline = values[baseField]
      ? calculateDeadline(values[baseField], daysToAdd)
      : "";
    setFieldValue(field, deadline);
  }, [baseField, daysToAdd, field, setFieldValue, values]);
};

const calculateDeadline = (startdate, daysToAdd) => {
  const deadline = new Date(startdate);
  deadline.setDate(deadline.getDate() + daysToAdd);

  /**
   * isoString is in the format of yyyy-mm-ddThh:mm:ss.sssZ
   * split string by 'T' and get the first element, which is the date
   */
  return deadline.toISOString().split("T")[0];
};
