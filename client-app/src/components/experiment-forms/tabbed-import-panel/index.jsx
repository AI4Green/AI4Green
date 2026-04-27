import {
  Box,
  FormLabel,
  TabList,
  TabPanel,
  TabPanels,
  Tabs,
  Text,
  useToast,
  VStack,
} from "@chakra-ui/react";
import { GLOBAL_PARAMETERS } from "constants";
import { useField } from "formik";
import { useCallback, useState } from "react";

import { FieldValueImporter } from "./field-value-importer";
import { RemovableTab } from "./removable-tab";
import { useParams } from "react-router-dom";

/**
 * Formik field for multiple tabs with import functionality.
 * Props:
 * - name: Name of the field in the formik form.
 * - isDisabled: Boolean to disable the field.
 * - sourceType: Type of the source. E.g. Plan or Lab Note
 * - fieldName: Name of the field to extract data from the source.
 * - label: Label for the field.
 * - buttonText: Text for the import button.
 * - Component: Component to render in the tab panel.
 * - isComponentDisabled: Boolean to disable the component.
 *
 * expects the field value to be an array of objects with the following structure:
 * - source: an object with properties.
 *  - type: Type of the source. E.g. Plan or Lab Note
 *  - name: Name/title if available. E.g. plan name or reaction name
 *  - id: Id of the source. E.g. Plan or Note id
 * - data: e.g. an array of product yield data.
 */
export const TabbedImportPanel = ({
  name,
  isDisabled,
  sourceType,
  fieldName,
  label,
  buttonText,
  Component,
  isComponentDisabled,
}) => {
  const { projectId } = useParams();
  const [field, , helpers] = useField(name);
  const toast = useToast();
  const [tabIndex, setTabIndex] = useState(0);

  const append = useCallback(
    (value) => {
      if (isSourceAlreadyAdded(value, field.value)) {
        toast({
          title: "Source already added",
          status: "warning",
          duration: GLOBAL_PARAMETERS.ToastDuration,
          isClosable: true,
          position: "top",
        });
        return;
      }

      setTabIndex(field.value.length);
      helpers.setValue([...field.value, value]);
    },
    [field.value, helpers, toast],
  );

  const removeTab = (index) => {
    helpers.setValue(field.value.filter((_, i) => i !== index));
    setTabIndex((currentTabIndex) =>
      currentTabIndex > 0 ? currentTabIndex - 1 : 0,
    );
  };

  return (
    <VStack
      spacing={4}
      w="full"
      p={4}
      border="1px"
      borderRadius="12"
      borderColor="gray.200"
    >
      <Box w="full">
        <FormLabel>{label}</FormLabel>
        <FieldValueImporter
          isDisabled={isDisabled}
          projectId={projectId}
          sourceType={sourceType}
          fieldName={fieldName}
          buttonText={buttonText}
          append={append}
        />
      </Box>

      {field?.value?.length > 0 ? (
        <Tabs
          index={tabIndex}
          onChange={setTabIndex}
          isLazy
          marginBottom={10}
          w="full"
        >
          <TabList>
            {field.value.map((val, index) => (
              <RemovableTab
                key={`${val.source?.type}.${val.source?.id}`}
                onRemove={() => removeTab(index)}
                variant="enclosed"
                isRemoveDisabled={isDisabled}
              >
                {val.source?.name}
              </RemovableTab>
            ))}
          </TabList>
          <TabPanels>
            {field.value.map((val, index) => (
              <TabPanel key={`${val.source?.type}.${val.source?.id}`}>
                <Component
                  name={`${name}[${index}].data`}
                  label={val.source?.name}
                  isDisabled={isDisabled || !isComponentDisabled}
                />
              </TabPanel>
            ))}
          </TabPanels>
        </Tabs>
      ) : (
        <Text pt={5} pb={5}>
          No data available
        </Text>
      )}
    </VStack>
  );
};

const isSourceAlreadyAdded = (value, existingValues) =>
  existingValues.some(
    (ex) =>
      ex.source.id === value.source.id && ex.source.type === value.source.type,
  );
