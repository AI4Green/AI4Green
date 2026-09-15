import { HStack, Icon, Text } from "@chakra-ui/react";
import { useTranslation } from "react-i18next";
import { FaExclamationTriangle } from "react-icons/fa";
import { string } from "yup";

import { FormikInput } from "./formik-input";

export const validationSchema = (t) => ({
  email: string()
    .email(t("validation.email_valid"))
    .required(t("validation.email_required")),
});

// Validate email against registration rules
export const validationSchemaRegRules = ({ t, validate }) => {
  return {
    email: validationSchema(t).email.test(
      "Valid email",
      "Invalid email or email not allowed for registration.",
      async (value) => {
        if (!value) return;
        const email = value.toLowerCase();
        const res = await validate(email);
        const { isValid } = await res.json();
        return isValid;
      },
    ),
  };
};

// Validate email based against registration rules and check if email is same as old one
export const validationSchemaExistingEmail = ({
  t,
  existingEmail,
  validate,
}) => ({
  email: validationSchemaRegRules({
    t,
    validate,
  }).email.notOneOf(
    [existingEmail],
    "Email must be different than the current one",
  ),
});

export const EmailField = ({ name = "email", hasCheckReminder, ...p }) => {
  const { t } = useTranslation();

  const checkReminder = (
    <HStack>
      <Icon as={FaExclamationTriangle} />
      <Text fontSize="xs">{t("fields.email_checkreminder")}</Text>
    </HStack>
  );

  return (
    <FormikInput
      name={name}
      type="email"
      isRequired
      label={t("fields.email")}
      placeholder={t("fields.email_placeholder")}
      fieldHelp={hasCheckReminder ? checkReminder : undefined}
      {...p}
    />
  );
};
