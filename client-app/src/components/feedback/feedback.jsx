import {
  HStack,
  Tag,
  TagCloseButton,
  TagLabel,
  TagLeftIcon,
  useDisclosure,
  useToast,
} from "@chakra-ui/react";
import { ActionButton } from "components/core/action-button";
import { LoadingIndicator } from "components/core/loading-indicator";
import {
  GLOBAL_PARAMETERS,
  SECTION_TYPES,
  SITE_ROLES,
  STAGES,
  STAGES_PERMISSIONS,
} from "constants";
import { useBackendApi, useUser } from "contexts";
import { useState } from "react";
import { useTranslation } from "react-i18next";
import { FaCheck, FaRegCommentAlt } from "react-icons/fa";

import { Comment } from "./comment";
import { CreateOrEditCommentModal } from "./comment/modal";

export const Feedback = ({ field, item }) => {
  const { user } = useUser();
  const { comments: action } = useBackendApi();

  const [isLoading, setIsLoading] = useState(false);

  const toast = useToast();
  const { t } = useTranslation();
  const { isOpen, onOpen, onClose } = useDisclosure();

  const isInstructor = user?.roles?.includes(SITE_ROLES.Instructor);

  if (!field || !item) return null;

  const handleApproval = async (isApproved) => {
    try {
      setIsLoading(true);
      const response = await action.setApprovalStatus(
        field.response.id,
        isApproved,
      );

      if (response && (response.status === 204 || response.status === 200)) {
        toast(
          toastOptions(
            isApproved ? "Approved" : "Approval cancelled",
            isApproved ? "success" : "warning",
          ),
        );
        await item.action.mutate();
      }
    } catch {
      toast(toastOptions(t("feedback.error_title"), "error"));
    } finally {
      setIsLoading(false);
    }
  };

  const canComment =
    isInstructor &&
    item.stage?.permissions.includes(STAGES_PERMISSIONS.InstructorCanComment);

  const canViewComment =
    (isInstructor || item.isOwner) &&
    item.stage?.name !== STAGES.Draft &&
    item.type !== SECTION_TYPES.ProjectGroup &&
    item.type !== SECTION_TYPES.Note &&
    field.feedback?.comments.total >= 1;

  const actions = {
    approve: {
      isEligible: () => canComment,
      label: "Approve",
      onClick: async () => await handleApproval(true),
      icon: <FaCheck />,
    },
    comment: {
      isEligible: () => canComment,
      label: "Add comment",
      onClick: onOpen,
      icon: <FaRegCommentAlt />,
    },
  };

  return (
    <HStack align="flex-start">
      {isLoading ? (
        <LoadingIndicator />
      ) : field.feedback?.approved ? (
        <Tag colorScheme="green" borderRadius="full" variant="outline">
          <TagLeftIcon as={FaCheck} />
          <TagLabel>Approved</TagLabel>
          {canComment && (
            <TagCloseButton onClick={() => handleApproval(false)} />
          )}
        </Tag>
      ) : (
        canComment && (
          <>
            <ActionButton actions={actions} size="xs" variant="outline" />
            {isOpen && (
              <CreateOrEditCommentModal
                fieldResponseId={field.response.id}
                isModalOpen={isOpen}
                onModalClose={onClose}
                mutate={item.action.mutate}
              />
            )}
          </>
        )
      )}
      {canViewComment && (
        <Comment
          field={field}
          canComment={canComment}
          mutate={item.action.mutate}
        />
      )}
    </HStack>
  );
};

export const toastOptions = (title, status) => ({
  position: "top",
  title,
  status,
  duration: GLOBAL_PARAMETERS.ToastDuration,
  isClosable: true,
});
