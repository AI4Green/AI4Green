import {
  Button,
  Modal as ChakraModal,
  ModalBody,
  ModalContent,
  ModalFooter,
  ModalHeader,
  ModalOverlay,
} from "@chakra-ui/react";
import { useCallback, useState } from "react";
import { FaRegCheckCircle, FaTimes } from "react-icons/fa";

export const Modal = ({
  size,
  title,
  body,
  isOpen,
  onClose,
  onAction,
  isLoading,
  actionBtnCaption = "Ok",
  actionBtnColorScheme = "green",
  actionBtnLeftIcon = <FaRegCheckCircle />,
  cancelBtnEnable = true,
  cancelBtnCaption = "Cancel",
  cancelBtnAction = onClose,
  closeOnOverlayClick = true,
  contentMaxW,
  contentMaxH,
  bodyMaxH,
  bodyOverflowY,
}) => (
  <ChakraModal
    closeOnEsc={closeOnOverlayClick}
    closeOnOverlayClick={closeOnOverlayClick}
    isOpen={isOpen}
    onClose={cancelBtnAction}
    size={size}
    isCentered
  >
    <ModalOverlay />
    <ModalContent maxW={contentMaxW} maxH={contentMaxH}>
      <ModalHeader fontSize="md" fontWeight="semibold">
        <Button
          onClick={cancelBtnAction}
          leftIcon={<FaTimes />}
          variant="ghost"
          size="sm"
          float="right"
        />

        {title}
      </ModalHeader>
      <ModalBody maxH={bodyMaxH} overflowY={bodyOverflowY}>
        {body}
      </ModalBody>
      <ModalFooter>
        {cancelBtnEnable && (
          <Button size="sm" onClick={cancelBtnAction} leftIcon={<FaTimes />}>
            {cancelBtnCaption}
          </Button>
        )}
        <Button
          size="sm"
          leftIcon={actionBtnLeftIcon}
          colorScheme={actionBtnColorScheme}
          onClick={onAction}
          ml={3}
          isLoading={isLoading}
        >
          {actionBtnCaption}
        </Button>
      </ModalFooter>
    </ModalContent>
  </ChakraModal>
);

export const useModalState = (location, navigate, formRef = null) => {
  const [isModalOpen, setIsModalOpen] = useState(false);
  const [isLoading, setIsLoading] = useState(false);
  const [feedback, setFeedback] = useState();

  const handleReset = useCallback(() => {
    setFeedback();
    setIsLoading(false);
    setIsModalOpen(false);
    formRef?.current?.resetForm();
    navigate(location.pathname, { replace: true });
  }, [location.pathname, navigate, formRef]);

  return {
    isModalOpen,
    setIsModalOpen,
    isLoading,
    setIsLoading,
    feedback,
    setFeedback,
    handleReset,
  };
};
