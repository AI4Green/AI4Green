import {
  HStack,
  Icon,
  Link,
  ListItem,
  Popover,
  PopoverBody,
  PopoverContent,
  PopoverTrigger,
  Text,
  UnorderedList,
} from "@chakra-ui/react";
import { useTranslation } from "react-i18next";
import { FaInfoCircle } from "react-icons/fa";
import { string } from "yup";

import { FormikInput } from "./formik-input";

export const minLength = 6;
export const validationSchema = (t) => ({
  password: string()
    .min(minLength, t("validation.password_minLength", { minLength }))
    .matches(/\d/, t("validation.password_digit"))
    .matches(/[A-Z]/, t("validation.password_uppercase"))
    .matches(/[a-z]/, t("validation.password_lowercase"))
    .matches(/[^A-Za-z0-9]/, t("validation.password_special"))
    .required(t("validation.password_required")),
});

const PasswordRequirementsTip = ({ minLength }) => {
  const { t } = useTranslation();

  const requirements = [
    t("validation.password_minLength", { minLength }),
    t("validation.password_digit"),
    t("validation.password_uppercase"),
    t("validation.password_lowercase"),
    t("validation.password_special"),
  ];

  return (
    <Popover returnFocusOnClose={false} usePortal>
      <PopoverTrigger>
        <Link fontSize="sm">
          <HStack align="center" spacing={1} mt={1}>
            <Icon as={FaInfoCircle} />
            <Text>{t("register.links.passwordRequirements")}</Text>
          </HStack>
        </Link>
      </PopoverTrigger>
      <PopoverContent bg="gray.50" borderColor="gray.200">
        <PopoverBody pl={4}>
          <UnorderedList>
            {requirements.map((x, i) => (
              <ListItem fontSize="xs" key={i}>
                {x}
              </ListItem>
            ))}
          </UnorderedList>
        </PopoverBody>
      </PopoverContent>
    </Popover>
  );
};

export const PasswordField = ({ name = "password", ...p }) => {
  const { t } = useTranslation();

  return (
    <FormikInput
      name={name}
      isRequired
      type="password"
      label={t("fields.password")}
      placeholder={t("fields.password")}
      fieldTip={<PasswordRequirementsTip minLength={minLength} />}
      {...p}
    />
  );
};
