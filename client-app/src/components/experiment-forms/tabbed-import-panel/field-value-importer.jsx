import {
  Alert,
  AlertIcon,
  Button,
  HStack,
  Icon,
  useDisclosure,
  VStack,
} from "@chakra-ui/react";
import { FormikInput, MultiSelectField } from "components/core/forms";
import { Modal } from "components/core/modal";
import { SECTION_TYPES, STAGES } from "constants";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { useEffect, useState } from "react";
import { useTranslation } from "react-i18next";
import { FaFileImport } from "react-icons/fa";
import { array, number, object } from "yup";

/**
 * This component is used to import field values from other sources (only Notes for now).
 * Renders a button that opens a modal to select a source and import the field values.
 * As modal is opened, it fetches the sources and the field, enabling the user to select a source and import the field response.
 * Fetched field id is used along with the source id to fetch the field response for the selected source.
 * The imported data is then passed to the append function.
 * Props:
 * - icon: Icon to display on the button.
 * - buttonText: Text for the button.
 * - isDisabled: Boolean to hide the button.
 * - projectId: Id of the project. Used to fetch the sources and field.
 * - sourceType: Type of the source. E.g. Lab Note. Used to fetch the sources list.
 * - fieldName: Name of the field to extract data from the source. Used to fetch the field.
 * - append: Function to append the imported data. To be configured in the parent component.
 */
export const FieldValueImporter = ({
  icon = <FaFileImport />,
  buttonText = "Import",
  isDisabled,
  projectId,
  sourceType,
  fieldName,
  append,
}) => {
  const { isOpen, onOpen, onClose } = useDisclosure();
  const modalProp = {
    isModalOpen: isOpen,
    onModalClose: onClose,
    buttonText: buttonText || "Import",
    projectId,
    sourceType,
    fieldName,
    append,
  };

  if (!isDisabled) {
    return (
      <>
        <Button
          isDisabled={isDisabled}
          onClick={onOpen}
          leftIcon={icon}
          size="sm"
          px={4}
          variant="outline"
        >
          {buttonText}
        </Button>
        <ImportModal {...modalProp} />
      </>
    );
  }
  return null;
};

// Modal component that renders the modal content.
const ImportModal = ({
  isModalOpen,
  onModalClose,
  buttonText,
  projectId,
  sourceType,
  fieldName,
  append,
}) => {
  const [isLoading, setIsLoading] = useState(false);
  const { action } = useSourceAction(sourceType);
  const { fields: fieldAction } = useBackendApi();

  const { sources, error: sourcesFetchError } = useFetchSources(
    isModalOpen,
    setIsLoading,
    projectId,
    sourceType,
    action.getList,
  );

  const { field, error: fieldFetchError } = useFetchField(
    isModalOpen,
    setIsLoading,
    projectId,
    sourceType,
    fieldName,
    fieldAction.getFieldByName,
  );

  const [feedback, setFeedback] = useState(
    sourcesFetchError || fieldFetchError,
  );

  const { t } = useTranslation();

  const handleSubmit = async (values) => {
    const sourceId = values.sourceId[0]; // only one selection allowed
    const sourceName = sources.find((source) => source.id === sourceId)?.name;
    try {
      setIsLoading(true);
      const response = await action.getFieldResponse(sourceId, field.id);
      const data = response.value;
      const source = {
        type: sourceType,
        id: sourceId,
        name: sourceName,
      };
      append({ source, data });
      await onModalClose();
    } catch (e) {
      console.error(e);
      switch (e?.response?.status) {
        case 404:
          setFeedback({
            status: "error",
            message: "Import failed. Data unavailable.",
          });
          break;
        default:
          setFeedback({
            status: "error",
            message: t("feedback.error_title"),
          });
          break;
      }
    } finally {
      setIsLoading(false);
    }
  };

  const modalBody = (
    <Form noValidate>
      <VStack align="stretch" spacing={4}>
        {feedback && (
          <Alert status={feedback.status}>
            <AlertIcon />
            {feedback.message}
          </Alert>
        )}
        <HStack spacing={4}>
          <Icon as={FaFileImport} color="blue.500" fontSize="5xl" />
          <VStack w="full">
            <FormikInput name="sourceType" label="Source Type" isDisabled />

            <MultiSelectField
              label="Source"
              placeholder="Select an item"
              name="sourceId"
              options={sources.map((source) => ({
                label: source.name,
                value: source.id,
              }))}
            />
          </VStack>
        </HStack>
      </VStack>
    </Form>
  );

  return (
    <Formik
      enableReinitialize
      initialValues={{
        sourceType,
        sourceId: [],
      }}
      onSubmit={handleSubmit}
      validationSchema={validationSchema(sources)}
    >
      {({ submitForm }) => (
        <Modal
          body={modalBody}
          title="Import"
          actionBtnCaption={buttonText}
          onAction={submitForm}
          actionBtnColorScheme="blue"
          isLoading={isLoading}
          isOpen={isModalOpen}
          onClose={
            feedback
              ? () => {
                  setFeedback(null);
                  onModalClose();
                }
              : onModalClose
          }
        />
      )}
    </Formik>
  );
};

const useSourceAction = (source) => {
  const { notes: noteAction } = useBackendApi();

  const sourceMap = {
    [SECTION_TYPES.Note]: {
      action: {
        getList: noteAction.getNotesList,
        getFieldResponse: noteAction.getNoteFieldResponse,
      },
    },
  };

  return sourceMap[source];
};

// Hook that fetches the sources list. e.g. Plans or Notes
const useFetchSources = (
  isModalOpen,
  setIsLoading,
  projectId,
  sourceType,
  getSources,
) => {
  const [sources, setSources] = useState([]);
  const [error, setError] = useState(null);

  useEffect(() => {
    if (isModalOpen) {
      setIsLoading(true);
      const fetchSources = async () => {
        try {
          const result = await getSources(projectId);
          const approvedSources = result
            .filter(
              (item) =>
                item.stage === STAGES.Approved ||
                item.plan.stage === STAGES.Approved,
            )
            .map((item) => ({
              id: item.id,
              name:
                item.title ||
                item.reactionName || // try to get the reaction name (note)
                `${sourceType} ${item.id}`,
            }));

          setSources(approvedSources);
          setIsLoading(false);
        } catch {
          setError("Failed to fetch sources");
          setIsLoading(false);
        }
      };
      fetchSources();
    }
  }, [getSources, isModalOpen, projectId, setIsLoading, sourceType]);

  return {
    sources,
    error: error && { status: "error", message: error },
  };
};

// Hook that fetches the field based on the field name.
const useFetchField = (
  isModalOpen,
  setIsLoading,
  projectId,
  sourceType,
  fieldName,
  getField,
) => {
  const [field, setField] = useState(null);
  const [error, setError] = useState(null);

  useEffect(() => {
    if (isModalOpen) {
      setIsLoading(true);
      const fetchField = async () => {
        try {
          const result = await getField(projectId, sourceType, fieldName);
          setField(result);
          setIsLoading(false);
        } catch {
          setError("Failed to fetch field");
          setIsLoading(false);
        }
      };
      fetchField();
    }
  }, [fieldName, getField, isModalOpen, projectId, setIsLoading, sourceType]);

  return {
    field,
    error: error && {
      status: "error",
      message: error,
    },
  };
};

const validationSchema = (sources) =>
  object().shape({
    sourceId: array()
      .min(1, "Please select a source")
      .of(
        number().oneOf(
          sources.map((source) => source.id),
          "Invalid source",
        ),
      )
      .required("Valid source required"),
  });
