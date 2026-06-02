import { Button, HStack, Text, useToast, VStack } from "@chakra-ui/react";
import { Breadcrumbs } from "components/core/breadcrumbs";
import { SectionField, validationSchema } from "components/section-field";
import { SectionHeader } from "components/section-header/header";
import {
  GLOBAL_PARAMETERS,
  SECTION_TYPES,
  STAGES_PERMISSIONS,
} from "constants";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { DefaultContentLayout } from "layouts/default";
import { useRef, useState } from "react";
import { useTranslation } from "react-i18next";
import { FaSave } from "react-icons/fa";

import { initialValues, prepareSubmissionData } from ".";

export const SectionForm = ({
  item,
  form,
  headerItems,
  breadcrumbItems,
  isInstructor,
}) => {
  const [isLoading, setIsLoading] = useState();
  const { t } = useTranslation();
  const toast = useToast();
  const formRef = useRef();
  const { sections: sectionApis } = useBackendApi();

  const sectionFields = form?.fields.map((x) => ({
    ...x,
    response: {
      id: x.fieldResponse?.[0]?.id || null,
      value: x.fieldResponse?.[0]?.fieldResponseValues?.at(-1)?.value || null,
    },
    feedback: x.fieldResponse?.[0]?.feedback || null,
  }));

  const handleSubmit = async (values, fields) => {
    const data = prepareSubmissionData(fields, values);

    const payload = {
      fieldResponses: JSON.stringify(data.fieldResponses), // used for updating existing field responses.
      newFieldResponses: JSON.stringify(data.newFieldResponses), // used for creating new field responses.
      sectionId: form.id,
      recordId: item.id,

      // file related data
      files: data.files,
      fileFieldResponses: JSON.stringify(data.fileFieldResponses),

      newFiles: data.newFiles,
      newFileFieldResponses: JSON.stringify(data.newFileFieldResponses),
    };

    // convert the payload to FormData
    const formData = new FormData();

    Object.entries(payload).forEach(([k, v]) => {
      if (k === "files" || k === "newFiles") {
        v.forEach((file) => formData.append(k, file)); // append files to the form
      } else if (Array.isArray(v)) {
        v.forEach((item) => formData.append(`${k}[]`, JSON.stringify(item)));
      } else if (typeof v === "object" && v !== null) {
        formData.append(k, JSON.stringify(v)); // stringify objects and append to the form
      } else {
        formData.append(k, v); // append the rest
      }
    });

    try {
      setIsLoading(true);
      const response = await item.action.save(formData);
      if (response && (response.status === 204 || response.status === 200)) {
        toast(toastOptions("Section values saved successfully", "success"));
        await item.action.mutate();
      }
    } catch (e) {
      console.log(e);
      toast(toastOptions(t("feedback.error_title"), "error"));
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <DefaultContentLayout>
      <Breadcrumbs items={breadcrumbItems} />
      <SectionHeader
        header={headerItems}
        project={{ name: form.reactionID }}
        action={
          <SectionFormAction
            item={item}
            isInstructor={isInstructor}
            isLoading={isLoading}
            formRef={formRef}
          />
        }
      />
      <Formik
        enableReinitialize
        initialValues={initialValues(
          sectionFields,
          item.id,
          form.id,
          sectionApis,
        )}
        validationSchema={validationSchema(sectionFields)}
        innerRef={formRef}
        onSubmit={async (values) => await handleSubmit(values, sectionFields)}
      >
        {({ values }) => (
          <Form noValidate>
            <VStack align="stretch" spacing={[3, 4]}>
              {sectionFields
                .sort((a, b) => a.sortOrder - b.sortOrder)
                .map(
                  (field) =>
                    !field.hidden && (
                      <SectionField
                        key={field.id}
                        fieldValues={values} // values is an collection of formik values, which can be accessed by using the field.id as key
                        field={field}
                        item={item}
                        fields={sectionFields}
                        isInstructor={isInstructor}
                        feedback={field.feedback}
                      />
                    ),
                )}
            </VStack>
          </Form>
        )}
      </Formik>
    </DefaultContentLayout>
  );
};

const toastOptions = (title, status) => ({
  position: "top",
  title,
  status,
  duration: GLOBAL_PARAMETERS.ToastDuration,
  isClosable: true,
});

const SectionFormAction = ({ item, isInstructor, isLoading, formRef }) => {
  console.log("Button Clicked");
  const hasRequiredPermissions = [
    STAGES_PERMISSIONS.OwnerCanEdit,
    STAGES_PERMISSIONS.OwnerCanEditCommented,
  ].some((permission) => item.approvalStatus?.permissions.includes(permission));

  const canUserSave =
    !isInstructor &&
    (item.type.toUpperCase() === SECTION_TYPES.ProjectGroup.toUpperCase() ||
      (item.isOwner && hasRequiredPermissions));

  return (
    <HStack pb={1}>
      {canUserSave && (
        <Button
          colorScheme="green"
          leftIcon={<FaSave />}
          size="sm"
          isLoading={isLoading}
          onClick={() => formRef.current.handleSubmit()}
        >
          <Text fontSize="sm" fontWeight="medium">
            Save
          </Text>
        </Button>
      )}
    </HStack>
  );
};
