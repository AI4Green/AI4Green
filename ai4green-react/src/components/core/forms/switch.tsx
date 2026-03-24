import {
  FormControl,
  HStack,
  Switch as ChakraSwitch,
  Text,
} from "@chakra-ui/react";
import { useField } from "formik";

export const Switch = ({ name, label, isRequired, ...p }) => {
  const [field, meta, helpers] = useField(name);
  return (
    <FormControl
      id={name}
      isRequired={isRequired}
      isInvalid={meta.error && meta.touched}
    >
      <HStack spacing={4} align="center">
        {label && <Text fontSize="sm">{label}</Text>}
        <ChakraSwitch
          size="sm"
          colorScheme="green"
          name={name}
          isChecked={field.value}
          onChange={(e) => helpers.setValue(e.target.checked)}
          {...p}
        />
      </HStack>
    </FormControl>
  );
};
