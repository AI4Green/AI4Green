import { FormErrorMessage, FormHelperText, Text } from "@chakra-ui/react";

/**
 * A reusable configurable component for displaying
 * help text and / or error text for a form control.
 *
 * MUST be inside a Chakra `<FormControl />`
 */
export const FormHelpError = ({
  isInvalid,
  help,
  collapseEmpty,
  replaceHelpWithError,
  error,
}) => {
  // collapseEmpty defaults to false so we
  // always display text containers (even if empty)
  // so the layout has a placeholder in case of errors,
  // preventing vertical shift when text appears.

  // help component, if any
  const displayHelp =
    collapseEmpty && !help ? null : (
      <FormHelperText lineHeight="normal">{help || <>&nbsp;</>}</FormHelperText>
    );

  // error component, if any
  const displayError =
    collapseEmpty && !error ? null : (
      <FormErrorMessage lineHeight="normal">
        <Text fontSize="xs">{error || <>&nbsp;</>}</Text>
      </FormErrorMessage>
    );

  // these are our actual rendered slots
  let helpSlot = replaceHelpWithError && isInvalid ? null : displayHelp;
  let errorSlot = replaceHelpWithError && !isInvalid ? null : displayError;

  return (
    <>
      {helpSlot}
      {errorSlot}
    </>
  );
};
