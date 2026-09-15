import {
  Alert,
  AlertIcon,
  HStack,
  Icon,
  Text,
  useToast,
  VStack,
} from "@chakra-ui/react";
import { Modal } from "components/core/modal";
import {
  GLOBAL_PARAMETERS,
  SECTION_TYPES,
  STAGE_TYPES,
  TITLE_ICON_COMPONENTS,
} from "constants";
import { useBackendApi } from "contexts";
import { useState } from "react";
import { useTranslation } from "react-i18next";

export const MoveStageModal = ({
  fixedNextStage,
  modalTitle,
  modalMessage,
  successMessage,
  isModalOpen,
  onModalClose,
  record,
  type,
  mutate,
}) => {
  const [isLoading, setIsLoading] = useState();
  const [feedback, setFeedback] = useState();

  const { t } = useTranslation();
  const toast = useToast();

  const { projectTypes, literatureReviews, plans, notes, reports } =
    useBackendApi();
  const { action, icon } = getStageItems(type, {
    projectTypes,
    literatureReviews,
    plans,
    notes,
    reports,
  });

  const handleAdvanceStage = async () => {
    try {
      setIsLoading(true);
      const response = await action.advanceStage(record.id, fixedNextStage);
      setIsLoading(false);

      if (response && (response.status === 204 || response.status === 200)) {
        toast({
          position: "top",
          title: successMessage,
          status: "success",
          duration: GLOBAL_PARAMETERS.ToastDuration,
          isClosable: true,
        });
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
      <HStack spacing={5}>
        <Icon as={icon} color={"green.500"} fontSize="5xl" />
        <VStack align="flex-end" flex={1}>
          <Text align="right">{modalMessage}</Text>
          <Text fontWeight="bold">{record.title}</Text>
        </VStack>
      </HStack>
    </VStack>
  );

  const resetState = () => {
    setFeedback();
    setIsLoading(false);
  };

  return (
    <Modal
      body={modalBody}
      title={modalTitle}
      actionBtnColorScheme="green"
      onAction={handleAdvanceStage}
      isLoading={isLoading}
      isOpen={isModalOpen}
      onClose={() => {
        resetState();
        onModalClose();
      }}
    />
  );
};

const getStageItems = (sectionType, apis) => {
  const { projectTypes, literatureReviews, plans, notes, reports } = apis;

  let items;
  switch (sectionType) {
    case STAGE_TYPES.ProjectType:
      items = {
        action: projectTypes,
        icon: TITLE_ICON_COMPONENTS.ProjectType,
      };
      break;
    case STAGE_TYPES.Plan:
      items = {
        action: plans,
        icon: TITLE_ICON_COMPONENTS.Plan,
      };
      break;
    case STAGE_TYPES.Report:
      items = {
        action: reports,
        icon: TITLE_ICON_COMPONENTS.Report,
      };
      break;
    case STAGE_TYPES.LiteratureReview:
      items = {
        action: literatureReviews,
        icon: TITLE_ICON_COMPONENTS.LiteratureReview,
      };
      break;
    case STAGE_TYPES.Note:
      items = {
        action: notes,
        icon: TITLE_ICON_COMPONENTS.Note,
      };
      break;
  }
  return items;
};
