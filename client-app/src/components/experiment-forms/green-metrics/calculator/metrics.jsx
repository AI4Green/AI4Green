import { Box, Button, FormLabel, HStack } from "@chakra-ui/react";
import { NumberInputField } from "components/core/forms";
import { Formik, useField } from "formik";
import { FaCalculator } from "react-icons/fa";

/**
 * This is a reusable calculator component used mainly for green metrics calculations
 * @param {string} name - field name (used for formik)
 * @param {string} title - calculator title
 * @param {Array} fields - array of fields object, such as name, label, isDisabled properties
 * @param {object} validationSchema - yup validation schema
 * @param {function} handleSubmit - function to handle form submission
 * @param {boolean} isDisabled - disable state
 */
export const MetricsCalculator = ({
  name,
  title,
  fields,
  validationSchema,
  handleSubmit,
  isDisabled,
}) => {
  const [field, , helpers] = useField(name);

  return (
    <Box w="full" align="flex-start">
      <FormLabel>{title}</FormLabel>
      <Formik
        enableReinitialize
        initialValues={{
          ...(field.value ||
            fields.reduce((acc, f) => ({ ...acc, [f.name]: 0 }), {})),
        }}
        onSubmit={(values) => handleSubmit(values, helpers)}
        validationSchema={validationSchema}
      >
        {({ handleSubmit }) => (
          <Box w="full" borderRadius={7} borderWidth={1} p={4}>
            <HStack align="flex-end">
              {fields.map((field) => (
                <NumberInputField
                  key={field.name}
                  name={field.name}
                  label={field.label}
                  placeholder={field.label}
                  isDisabled={field.isDisabled || isDisabled}
                  isRequired
                />
              ))}
            </HStack>
            <Button
              onClick={handleSubmit}
              colorScheme="green"
              size="sm"
              px={4}
              leftIcon={<FaCalculator />}
              hidden={isDisabled}
            >
              Calculate
            </Button>
          </Box>
        )}
      </Formik>
    </Box>
  );
};
