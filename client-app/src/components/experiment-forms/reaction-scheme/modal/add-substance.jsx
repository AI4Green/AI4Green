import {
  Box,
  FormLabel,
  HStack,
  Icon,
  Text,
  useToast,
  VStack,
} from "@chakra-ui/react";
// import { useSolventsList } from "api";
import { Modal } from "components/core/modal";
import { GLOBAL_PARAMETERS } from "constants";
import { useBackendApi } from "contexts";
import { ErrorMessage, Form, Formik } from "formik";
import { useRef } from "react";
import { FaFlask, FaVial } from "react-icons/fa";
import AsyncSelect from "react-select/async";
import { theme } from "themes/theme";
import { object, string } from "yup";

export const AddSubstanceModal = ({
  isModalOpen,
  onModalClose,
  setTableData,
  isAddingSolvent,
}) => {
  const formRef = useRef();
  const toast = useToast();

  const { reactionTable: action } = useBackendApi();
  const { data: solvents } = useSolventsList();

  const solventsOptions = solvents?.map((item) => ({
    value: item.name,
    label: item.name,
    isSolvent: true,
    flag: item.flag,
  }));

  const handleAddSubstance = async (values) => {
    try {
      const data = isAddingSolvent
        ? await action.getSolvent(values.substance)
        : await action.getReagent(values.substance);
      setTableData((old) => [...old, createRowData(data, values)]);
      onModalClose();
    } catch (error) {
      toast({
        title: "An error occurred.",
        description: error.message,
        status: "error",
        duration: GLOBAL_PARAMETERS.ToastDurationLong,
        isClosable: true,
        direction: "top",
      });
    }
  };

  const timeoutRef = useRef(null);

  const modalBody = (
    <Formik
      enableReinitialize
      innerRef={formRef}
      initialValues={{
        substanceType: !isAddingSolvent ? "Reagent" : "Solvent",
        substance: "",
      }}
      onSubmit={handleAddSubstance}
      validationSchema={validationSchema}
    >
      {({ setFieldValue }) => (
        <Form noValidate>
          <VStack align="stretch" spacing={4}>
            <HStack spacing={5}>
              <Icon
                as={isAddingSolvent ? FaVial : FaFlask}
                color={isAddingSolvent ? "teal" : "pink.600"}
                fontSize="5xl"
              />
              <Box flex={1}>
                <FormLabel>Substance</FormLabel>
                <AsyncSelect
                  cacheOptions
                  loadOptions={(inputValue, callback) =>
                    isAddingSolvent
                      ? callback(
                          solventsOptions.filter((option) =>
                            option.label
                              .toLowerCase()
                              .includes(inputValue.toLowerCase()),
                          ),
                        )
                      : loadCompounds(
                          inputValue,
                          callback,
                          timeoutRef,
                          action.getCompounds,
                        )
                  }
                  defaultOptions={isAddingSolvent ? solventsOptions : []}
                  placeholder="Start typing to search for a substance"
                  onChange={(option) => {
                    setFieldValue("substance", option?.value || "");
                  }}
                  styles={{
                    option: (provided, state) => ({
                      ...provided,
                      backgroundColor: state.data.flag
                        ? state.isFocused
                          ? provided.backgroundColor
                          : colorMap[state.data.flag].bgColor
                        : provided.backgroundColor,
                      color: state.data.flag
                        ? state.isFocused
                          ? provided.color
                          : colorMap[state.data.flag].color
                        : provided.color,
                      fontWeight: state.isFocused && theme.fontWeights.medium,
                    }),
                  }}
                />
                <ErrorMessage name="substance">
                  {(msg) => (
                    <Text fontSize="sm" color="red.500">
                      {msg}
                    </Text>
                  )}
                </ErrorMessage>
              </Box>
            </HStack>
          </VStack>
        </Form>
      )}
    </Formik>
  );
  return (
    <Modal
      body={modalBody}
      title={`Add ${isAddingSolvent ? "Solvent" : "Reagent"}`}
      actionBtnCaption="Add"
      onAction={() => formRef.current.handleSubmit()}
      isOpen={isModalOpen}
      onClose={onModalClose}
      actionBtnColorScheme={isAddingSolvent ? "teal" : "pink"}
    />
  );
};

// using timeout to debounce search instead of making request every key stroke
const loadCompounds = (inputValue, callback, timeoutRef, getCompounds) => {
  if (timeoutRef.current) clearTimeout(timeoutRef.current);

  timeoutRef.current = setTimeout(async () => {
    if (!inputValue) {
      callback([]);
      return;
    }

    try {
      const response = await getCompounds(inputValue);
      const options = response?.map((item) => ({
        value: item.name,
        label: item.name,
      }));
      callback(options);
    } catch (error) {
      console.error(error);
      callback([]);
    }
  }, 500);
};

const createRowData = (data, values) => ({
  manualEntry: true,
  substanceType: values?.substanceType,
  substancesUsed: data?.name,
  molWeight: data?.molecularWeight,
  density: data?.density,
  hazards: data?.hazards,
  mass: { value: 0, unit: "" },
});

const validationSchema = object().shape({
  substance: string().required("Please select a substance"),
  substanceType: string().required("Substance type is required"),
});

// color map extracted from ai4green
const colorMap = {
  1: { bgColor: "#8b0000", color: "#FFFFFF", rate: "hazard-highly-hazardous" },
  4: { bgColor: "#00ff00", color: "#000000", rate: "hazard-acceptable" },
  3: { bgColor: "#ffff00", color: "#000000", rate: "hazard-warning" },
  2: { bgColor: "#ff0000", color: "#000000", rate: "hazard-hazardous" },
  5: { bgColor: "#FFFFFF", color: "#000000", rate: "non-chem21" },
};
