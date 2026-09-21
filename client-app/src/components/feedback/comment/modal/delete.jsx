import { Alert, AlertIcon, Text, useToast, VStack } from "@chakra-ui/react";
import { Modal } from "components/core/modal";
import { toastOptions } from "components/feedback/feedback";
import { useBackendApi } from "contexts";
import { useState } from "react";
import { useTranslation } from "react-i18next";

export const DeleteCommentModal = ({
  isModalOpen,
  onModalClose,
  comment,
  mutate,
}) => {
  const { comments: action } = useBackendApi();

  const [isLoading, setIsLoading] = useState();
  const [feedback, setFeedback] = useState();

  const { t } = useTranslation();
  const toast = useToast();

  const handleDelete = async () => {
    try {
      setIsLoading(true);
      const response = await action.delete(comment.id);
      setIsLoading(false);

      if (response && (response.status === 204 || response.status === 200)) {
        toast(toastOptions("Comment deleted", "success"));
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

  const modalBody = (
    <VStack>
      {feedback && (
        <Alert status={feedback.status}>
          <AlertIcon />
          {feedback.message}
        </Alert>
      )}
      <Text>Are you sure you want to delete this comment?</Text>
      <Text fontWeight="bold">{comment.value}</Text>
    </VStack>
  );

  const resetState = () => {
    setFeedback();
    setIsLoading(false);
  };

  return (
    <Modal
      body={modalBody}
      title="Delete Confirmation"
      actionBtnCaption="Delete"
      actionBtnColorScheme="red"
      onAction={handleDelete}
      isLoading={isLoading}
      isOpen={isModalOpen}
      onClose={() => {
        resetState();
        onModalClose();
      }}
    />
  );
};
