import {
  FormControl,
  FormLabel,
  IconButton,
  Input,
  InputGroup,
  InputLeftElement,
  InputRightElement,
  useToast,
} from "@chakra-ui/react";
import { useField } from "formik";
import { useDebounce } from "helpers/hooks";
import { useEffect, useState } from "react";
import { FaEye, FaEyeSlash, FaRegCopy } from "react-icons/fa";

import { FormHelpError } from "./form-help-error";

export const FormikInput = ({
  name,
  label,
  placeholder,
  type = "text",
  isRequired,
  fieldTip,
  fieldHelp,
  collapseError,
  ...p
}) => {
  const toast = useToast();
  const [field, meta, helpers] = useField({ name, type });

  const [isMasked, setIsMasked] = useState(type === "password");

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

  const inputField = (
    <Input
      type={isMasked ? "password" : type === "password" ? "text" : type}
      placeholder={placeholder}
      {...p}
      value={value}
      onChange={handleChange}
      onBlur={type === "email" ? field.onBlur : undefined}
      _disabled={{ opacity: 0.7 }}
    />
  );

  const onClickCopyToClipboard = (value) => {
    // handle copy to clipboard action
    navigator.clipboard.writeText(value);
    toast({
      position: "top",
      title: "Copied to clipboard",
      status: "success",
      duration: 700,
      isClosable: true,
    });
  };

  return (
    <FormControl
      id={field.name}
      isRequired={isRequired}
      isInvalid={meta.error && meta.touched}
    >
      {label && <FormLabel>{label}</FormLabel>}

      {type === "password" || type === "readOnly" ? (
        <InputGroup>
          {inputField}
          {type == "password" ? (
            <InputLeftElement>
              <IconButton
                variant="ghost"
                onClick={() => setIsMasked(!isMasked)}
                size="sm"
                icon={isMasked ? <FaEye /> : <FaEyeSlash />}
              />
            </InputLeftElement>
          ) : (
            value && ( // display copy icon to allow user to copy value to clipboard
              <InputRightElement>
                <IconButton
                  variant="ghost"
                  onClick={() => onClickCopyToClipboard(value)}
                  size="sm"
                  icon={<FaRegCopy />}
                />
              </InputRightElement>
            )
          )}
        </InputGroup>
      ) : (
        inputField
      )}

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
