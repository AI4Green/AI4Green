import {
  Box,
  HStack,
  Popover,
  PopoverArrow,
  PopoverBody,
  PopoverCloseButton,
  PopoverContent,
  PopoverHeader,
  PopoverTrigger,
  Portal,
  Text,
  VStack,
} from "@chakra-ui/react";
import { LoadingIndicator } from "components/core/loading-indicator";
import { NotificationBadge } from "components/core/notification-badge";
import { SITE_ROLES } from "constants";
import { useBackendApi, useUser } from "contexts";
import { getFormattedDate } from "helpers";
import { useState } from "react";
import { FaRegCommentDots } from "react-icons/fa";

import { CommentActions } from "./actions";

export const Comment = ({ field, canComment, mutate }) => {
  const { comments: action } = useBackendApi();

  const [isLoading, setIsLoading] = useState(false);
  const [feedback, setFeedback] = useState();
  const [commentLogs, setCommentLogs] = useState([]);

  const handleOpen = async () => {
    try {
      setIsLoading(true);
      const data = await action.getCommentLogs(field.response.id);
      setCommentLogs(data);
      feedback && setFeedback(null);
    } catch (error) {
      console.error(error);
      setFeedback("Error fetching comments");
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <Popover onOpen={handleOpen}>
      <PopoverTrigger>
        <Box>
          <NotificationBadge
            as="button"
            type="button"
            count={field.feedback?.comments.unread}
            icon={FaRegCommentDots}
            aria-label="Comments"
          />
        </Box>
      </PopoverTrigger>
      <Portal>
        <PopoverContent>
          <PopoverArrow />
          <PopoverHeader>Comments</PopoverHeader>
          <PopoverCloseButton />
          <PopoverBody overflowY="auto" maxH="300px">
            {isLoading ? (
              <LoadingIndicator />
            ) : feedback ? (
              <Text color="red">{feedback}</Text>
            ) : (
              commentLogs
                ?.sort(
                  (a, b) => new Date(b.commentDate) - new Date(a.commentDate),
                )
                .map((log) => (
                  <CommentLog
                    key={log.id}
                    comment={log}
                    fieldResponseId={field.response.id}
                    canComment={canComment}
                    mutate={mutate}
                  />
                ))
            )}
          </PopoverBody>
        </PopoverContent>
      </Portal>
    </Popover>
  );
};

const CommentLog = ({ comment, fieldResponseId, canComment, mutate }) => {
  const { user } = useUser();
  const isInstructor = user?.roles?.includes(SITE_ROLES.Instructor);

  return (
    <VStack
      align="stretch"
      borderBottomWidth={1}
      borderRadius={4}
      p={1}
      my={2}
      fontSize="sm"
      bgColor={!comment.read && "gray.50"}
    >
      <Box display="flex" justifyContent="flex-end" mb={-3}>
        <CommentActions
          comment={comment}
          fieldResponseId={fieldResponseId}
          canComment={canComment}
          isInstructor={isInstructor}
          mutate={mutate}
        />
      </Box>

      <Text>{comment.value}</Text>
      <HStack fontSize="xxs" justify="flex-end" color="gray.600">
        <Text fontWeight="medium">{comment.owner}</Text>
        <Text>{getFormattedDate(comment.commentDate, user.uiCulture)}</Text>
      </HStack>
    </VStack>
  );
};
