import { useToast } from "@chakra-ui/react";
// import { useProject } from "api";
import { Modal, useModalState } from "components/core/modal";
import { InviteModal } from "components/project/modal";
import { GLOBAL_PARAMETERS, TITLE_ICON_COMPONENTS } from "constants";
import { useBackendApi } from "contexts";
import { useEffect, useRef } from "react";
import { useTranslation } from "react-i18next";
import {
  useLocation,
  useNavigate,
  useParams,
  useSearchParams,
} from "react-router-dom";

export const InstructorInviteModal = ({ mutate }) => {
  const [searchParams] = useSearchParams();
  const { projectId } = useParams();
  const { data: project } = useProject(projectId);
  const isInviteAction = searchParams.get("action") === "invite-instructors";

  const navigate = useNavigate();
  const location = useLocation();

  const { projects: action } = useBackendApi();

  const { t } = useTranslation();
  const toast = useToast();

  const formRef = useRef();

  const {
    isModalOpen,
    setIsModalOpen,
    isLoading,
    setIsLoading,
    feedback,
    setFeedback,
    handleReset,
  } = useModalState(location, navigate, formRef);

  useEffect(() => {
    setIsModalOpen(true);
  }, [projectId, isInviteAction, setIsModalOpen]);

  const handleSubmit = async ({ emails }) => {
    try {
      setIsLoading(true);
      const response = await action.inviteInstructors(projectId, {
        emails,
      });
      setIsLoading(false);

      if (response && (response.status === 204 || response.status === 200)) {
        toast({
          title: `Instructors invited successfully`,
          status: "success",
          duration: GLOBAL_PARAMETERS.ToastDuration,
          position: "top",
          isClosable: true,
        });
        handleReset();
        await mutate();
      }
    } catch {
      setFeedback({
        status: "error",
        message: t("feedback.error_title"),
      });
    }
  };

  if (!project) {
    navigate(location.pathname, { replace: true });
    return null;
  }

  return (
    <Modal
      body={
        <InviteModal
          ref={formRef}
          handleSubmit={handleSubmit}
          feedback={feedback}
          title="Invite instructors to the project."
          tags={[
            {
              label: project.name,
              colorScheme: "green",
              leftIcon: TITLE_ICON_COMPONENTS.Project,
            },
          ]}
        />
      }
      title="Project instructor invitation"
      actionBtnCaption="Invite"
      onAction={() => formRef.current.handleSubmit()}
      isLoading={isLoading}
      isOpen={isModalOpen}
      onClose={handleReset}
    />
  );
};
