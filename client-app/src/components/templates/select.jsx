import { Form, Formik } from "formik";
import { Alert, AlertIcon } from "@chakra-ui/react";
import { MultiSelectField } from "components/core/forms";

export const TemplateSelectorForm = ({
  formRef,
  templateLabel,
  templates,
  initialValues,
  feedback,
  onSubmit,
}) => {
  const templateOptions = templates.map((template) => ({
    value: template.id,
    label: template.name,
  }));

  return (
    <Formik
      innerRef={formRef}
      initialValues={initialValues}
      onSubmit={onSubmit}
    >
      <Form>
        {feedback && (
          <Alert status={feedback.status}>
            <AlertIcon />
            {feedback.message}
          </Alert>
        )}

        <MultiSelectField
          name="templateId"
          label={templateLabel}
          options={templateOptions}
        />
      </Form>
    </Formik>
  );
};
