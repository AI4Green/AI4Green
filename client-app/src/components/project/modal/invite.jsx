import {
  Alert,
  AlertIcon,
  HStack,
  Icon,
  Text,
  Textarea,
  Tooltip,
  VStack,
} from "@chakra-ui/react";
import { Badge } from "components/core/Badge";
import { MultiSelectField } from "components/core/forms";
import { Form, Formik } from "formik";
import { useState } from "react";
import { IoMdInformationCircleOutline } from "react-icons/io";
import { array, object, string } from "yup";

export const InviteModal = ({ ref, title, tags, handleSubmit, feedback }) => {
  const [emailList, setEmailList] = useState([]);

  return (
    <Formik
      enableReinitialize
      innerRef={ref}
      initialValues={{
        emails: [],
      }}
      onSubmit={handleSubmit}
      validationSchema={validationSchema}
    >
      {({ setFieldValue }) => {
        const handleEmailTextAreaChange = ({ target: { value } }) => {
          const emailSchema = string().email("Invalid email format");

          const list = value.split(",").map((item) => item.trim());
          const uniqueList = [...new Set(list)];
          const validEmailList = uniqueList.filter((item) => {
            return item !== "" && emailSchema.isValidSync(item);
          });

          setEmailList(validEmailList);
          setFieldValue("emails", validEmailList);
        };
        return (
          <Form noValidate>
            <VStack align="flex-start" spacing={8}>
              {feedback && (
                <Alert status={feedback.status}>
                  <AlertIcon />
                  {feedback.message}
                </Alert>
              )}

              <VStack align="flex-start" spacing={2}>
                {title && (
                  <Text fontSize="sm" fontWeight="light" color="gray.600">
                    {title}
                  </Text>
                )}
                {tags?.length > 0 && (
                  <HStack>
                    {tags.map((tag) => (
                      <Badge
                        key={tag.label}
                        colorScheme={tag.colorScheme}
                        label={tag.label}
                        leftIcon={tag.leftIcon}
                        variant="outline"
                        fontSize="xxs"
                      />
                    ))}
                  </HStack>
                )}
              </VStack>

              <HStack w="full">
                <Textarea
                  placeholder="Add emails for invitation"
                  rows="8"
                  onChange={(value) =>
                    handleEmailTextAreaChange(value, setFieldValue)
                  }
                />
                <Tooltip
                  label="Enter comma-separated email addresses. Invalid emails and duplicates will be automatically filtered out. Valid emails will appear in the selection field below where you can review and remove any before sending invitations."
                  placement="right"
                  hasArrow
                  fontSize="xs"
                >
                  <Icon
                    as={IoMdInformationCircleOutline}
                    boxSize={4}
                    cursor="help"
                    color="gray.500"
                  />
                </Tooltip>
              </HStack>

              <MultiSelectField
                isRequired
                isMulti
                label="Emails"
                placeholder="Select emails"
                name="emails"
                options={emailList.map((email) => ({
                  label: email,
                  value: email,
                }))}
              />
            </VStack>
          </Form>
        );
      }}
    </Formik>
  );
};

const validationSchema = () =>
  object().shape({
    emails: array()
      .required("Emails list required")
      .min(1, "Please select at least one email"),
  });
