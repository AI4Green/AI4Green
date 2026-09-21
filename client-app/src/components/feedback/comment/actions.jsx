import { IconButton, useDisclosure, useToast } from "@chakra-ui/react";
import { toastOptions } from "components/feedback/feedback";
import { useBackendApi } from "contexts";
import { FaEdit, FaRegDotCircle, FaRegTimesCircle } from "react-icons/fa";

import { CreateOrEditCommentModal, DeleteCommentModal } from "./modal";

export const CommentActions = ({
  comment,
  fieldResponseId,
  canComment,
  isInstructor,
  mutate,
}) => {
  const { comments: action } = useBackendApi();
  const toast = useToast();

  const editState = useDisclosure();
  const deleteState = useDisclosure();

  const handleMarkCommentAsRead = async () => {
    try {
      const response = await action.markAsRead(comment.id);
      if (response && response.status === 204) {
        toast(toastOptions("Comment marked as read", "success"));
        await mutate();
      }
    } catch {
      toast(toastOptions("Error", "error"));
    }
  };

  const showMarkAsRead = !comment.read && !isInstructor;

  return (
    <>
      {showMarkAsRead && (
        <IconButton
          icon={<FaRegDotCircle />}
          isRound
          size="xs"
          variant="ghost"
          onClick={handleMarkCommentAsRead}
        />
      )}
      {canComment && (
        <>
          <IconButton
            icon={<FaEdit />}
            isRound
            size="xs"
            variant="ghost"
            onClick={editState.onOpen}
          />
          {editState.isOpen && (
            <CreateOrEditCommentModal
              comment={comment}
              fieldResponseId={fieldResponseId}
              isModalOpen={editState.isOpen}
              onModalClose={editState.onClose}
              mutate={mutate}
            />
          )}

          <IconButton
            icon={<FaRegTimesCircle />}
            isRound
            size="xs"
            variant="ghost"
            onClick={deleteState.onOpen}
          />
          {deleteState.isOpen && (
            <DeleteCommentModal
              comment={comment}
              fieldResponseId={fieldResponseId}
              isModalOpen={deleteState.isOpen}
              onModalClose={deleteState.onClose}
              mutate={mutate}
            />
          )}
        </>
      )}
    </>
  );
};
