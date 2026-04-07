import {
  Badge,
  FormControl,
  FormLabel,
  Textarea,
  VStack,
} from "@chakra-ui/react";
import { useField } from "formik";
import { countWords } from "helpers";
import { useDebounce } from "helpers/hooks";
import { useEffect, useState } from "react";

import { FormHelpError } from "./form-help-error";

const WordCountBadge = ({ value, limit }) => {
  const count = countWords(value);
  return limit ? (
    <Badge colorScheme={count > limit ? "red" : undefined}>
      Word Count: {count} / {limit}
    </Badge>
  ) : null;
};

export const TextAreaField = ({
  name,
  title,
  wordLimit,
  placeholder,
  collapseError,
  fieldTip,
  fieldHelp,
  isRequired,
  ...p
}) => {
  const [field, meta, helpers] = useField({ name, type: "text" });

  const [value, setValue] = useState(field.value);
  const debouncedValue = useDebounce(value, 150);

  const handleChange = ({ target: { value } }) => {
    setValue(value);
  };

  useEffect(() => {
    helpers.setValue(debouncedValue);
  }, [debouncedValue, helpers]);

  useEffect(() => {
    setValue(field.value);
  }, [field.value]);

  return (
    <VStack w="100%" align="start">
      <FormControl
        id={field.name}
        isRequired={isRequired}
        isInvalid={meta.error && meta.touched}
      >
        <FormLabel>{title}</FormLabel>
        <Textarea
          value={value}
          placeholder={placeholder}
          onChange={handleChange}
          rows="6"
          _focus={{ borderColor: "brand.100" }}
          _disabled={{ opacity: 0.7 }}
          {...p}
        />

        {wordLimit && <WordCountBadge value={value} limit={wordLimit} />}

        {fieldTip}

        <FormHelpError
          isInvalid={meta.error && meta.touched}
          error={meta.error}
          help={fieldHelp}
          collapseEmpty={collapseError}
          replaceHelpWithError
        />
      </FormControl>
    </VStack>
  );
};
