import { Alert, AlertIcon, useToast, VStack } from "@chakra-ui/react";
import { TextAreaField } from "components/core/forms";
import { Modal } from "components/core/modal";
import { toastOptions } from "components/feedback/feedback";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { useRef, useState } from "react";
import { useTranslation } from "react-i18next";
import { object, string } from "yup";

export const CreateOrEditCommentModal = ({
  fieldResponseId,
  isModalOpen,
  onModalClose,
  comment,
  mutate,
}) => {
  const { comments: action } = useBackendApi();

  const [isLoading, setIsLoading] = useState();
  const [feedback, setFeedback] = useState();

  const toast = useToast();
  const formRef = useRef();
  const { t } = useTranslation();

  const handleSubmit = async (values) => {
    try {
      setIsLoading(true);
      const response = !comment
        ? await action.create(values)
        : await action.edit({ ...values, read: false }, comment?.id);
      setIsLoading(false);

      if (response && response.status === 200) {
        toast(
          toastOptions(
            `${comment?.id ? "Comment updated" : "Comment created"}`,
            "success",
          ),
        );

        await mutate();
        onModalClose();
      }
    } catch {
      setFeedback({
        status: "error",
        message: t("feedback.error_title"),
      });
    }
  };
  const validationSchema = () =>
    object().shape({
      value: string().required("Comment value is required"),
    });

  const modalBody = (
    <Formik
      enableReinitialize
      innerRef={formRef}
      initialValues={{
        fieldResponseId,
        value: comment?.value || "",
      }}
      onSubmit={handleSubmit}
      validationSchema={validationSchema()}
    >
      <Form noValidate>
        <VStack align="stretch" spacing={4}>
          {feedback && (
            <Alert status={feedback.status}>
              <AlertIcon />
              {feedback.message}
            </Alert>
          )}
          <TextAreaField
            name="value"
            title="Comment"
            placeholder="Enter your comment here"
            isRequired
          />
        </VStack>
      </Form>
    </Formik>
  );

  const resetState = () => {
    setFeedback();
    setIsLoading(false);
  };

  return (
    <Modal
      body={modalBody}
      title={`${comment?.id ? "Edit" : "Create"} Comment`}
      actionBtnCaption="Comment"
      onAction={() => formRef.current.handleSubmit()}
      isLoading={isLoading}
      isOpen={isModalOpen}
      onClose={() => {
        resetState();
        onModalClose();
      }}
    />
  );
};
