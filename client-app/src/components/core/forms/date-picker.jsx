import { FormControl, FormLabel, Input } from "@chakra-ui/react";
import { useField } from "formik";

import { FormHelpError } from "./form-help-error";

/**
 *
 * Note: This date value is in the format of yyyy-mm-dd
 * https://developer.mozilla.org/en-US/docs/Web/HTML/Element/input/date
 */

export const Datepicker = ({ name, label, isRequired }) => {
  const [field, meta, helpers] = useField(name);

  return (
    <FormControl
      isRequired={isRequired}
      id={field.name}
      isInvalid={meta.error && meta.touched}
    >
      <FormLabel>{label}</FormLabel>
      <Input
        placeholder="Select Date"
        size="md"
        type="date"
        value={field.value}
        onChange={(e) => helpers.setValue(e.target.value)}
      />

      <FormHelpError
        isInvalid={meta.touched && meta.error}
        error={meta.error}
        collapseEmpty
        replaceHelpWithError
      />
    </FormControl>
  );
};
