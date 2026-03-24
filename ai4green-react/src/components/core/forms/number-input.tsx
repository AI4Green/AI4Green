import {
  Box,
  FormControl,
  FormLabel,
  NumberDecrementStepper,
  NumberIncrementStepper,
  NumberInput,
  NumberInputField as ChakraNumberInputField,
  NumberInputStepper,
} from "@chakra-ui/react";
import { useField } from "formik";
import { useDebounce } from "helpers/hooks";
import { useEffect, useState } from "react";

import { FormHelpError } from "./form-help-error";

export const NumberInputField = ({
  name,
  label,
  placeholder,
  isRequired,
  fieldTip,
  fieldHelp,
  collapseError,
  isDisabled,
}) => {
  const [field, meta, helpers] = useField(name);
  const [value, setValue] = useState(field.value);
  const debouncedValue = useDebounce(value, 150);

  const handleChange = (value) => {
    setValue(parseFloat(value));
  };
  useEffect(() => {
    helpers.setValue(debouncedValue);
  }, [debouncedValue, helpers]);

  useEffect(() => {
    setValue(field.value);
  }, [field.value]);

  return (
    <FormControl
      id={field.name}
      isRequired={isRequired}
      isInvalid={meta.error && meta.touched}
    >
      {label && <FormLabel>{label}</FormLabel>}

      <Box maxW={48}>
        <NumberInput
          size="sm"
          step={0.2}
          onChange={handleChange}
          value={value}
          isDisabled={isDisabled}
        >
          <ChakraNumberInputField placeholder={placeholder} />
          <NumberInputStepper>
            <NumberIncrementStepper />
            <NumberDecrementStepper />
          </NumberInputStepper>
        </NumberInput>
      </Box>
      {fieldTip}

      <FormHelpError
        isInvalid={meta.error && meta.touched}
        error={meta.error}
        help={fieldHelp}
        collapseEmpty={collapseError}
        replaceHelpWithError
      />
    </FormControl>
  );
};
