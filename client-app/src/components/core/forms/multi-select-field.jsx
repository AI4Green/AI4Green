/*
  Select component, which can be used as either single select or multi select
*/

import {
  Box,
  FormControl,
  FormLabel,
  HStack,
  Icon,
  Text,
  theme,
  Tooltip,
} from "@chakra-ui/react";
import { Field } from "formik";
import { IoMdInformationCircleOutline } from "react-icons/io";
import Select from "react-select";

import { FormHelpError } from "./form-help-error";

export const MultiSelectField = ({
  options,
  name,
  placeholder,
  label,
  isMulti,
  initialValue,
  isDisabled,
  isRequired,
}) => (
  <Field name={name} initialValue={initialValue}>
    {({ field, form }) => {
      const onChange = (option) => {
        const selectedValue = isMulti
          ? option.map((item) => item.value)
          : [option?.value || ""];
        form.setFieldValue(field.name, selectedValue);
      };

      const selected = () => {
        return options
          ? options.filter((option) => field.value.indexOf(option.value) >= 0)
          : [];
      };

      return (
        <FormControl
          isRequired={isRequired}
          isInvalid={form.errors[field.name]}
        >
          <FormLabel>{label}</FormLabel>
          <Select
            name={field.name}
            value={selected()}
            onChange={onChange}
            placeholder={placeholder}
            options={options}
            isMulti={isMulti}
            isDisabled={isDisabled}
            components={{ Option: CustomOption }}
          />
          <FormHelpError
            isInvalid={form.errors[field.name] && form.touched[field.name]}
            error={form.errors[field.name]}
            collapseEmpty
          />
        </FormControl>
      );
    }}
  </Field>
);

const CustomOption = ({ data, ...props }) => {
  return (
    <Box {...props.innerProps} _hover={{ bg: theme.colors.gray[100] }}>
      <HStack
        fontSize="sm"
        fontWeight="thin"
        py={1}
        px={2}
        justify="space-between"
      >
        <Text>{data.label}</Text>
        {data.description && (
          <Tooltip
            label={data.description}
            placement="right"
            hasArrow
            fontSize="xs"
          >
            <Icon
              as={IoMdInformationCircleOutline}
              ml={2}
              boxSize={4}
              cursor="help"
            />
          </Tooltip>
        )}
      </HStack>
    </Box>
  );
};
