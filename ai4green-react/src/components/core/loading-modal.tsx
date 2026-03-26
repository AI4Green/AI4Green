import { Modal, ModalBody, ModalContent, ModalOverlay } from "@chakra-ui/react";

import { LoadingIndicator } from "./loading-indicator";

export const LoadingModal = ({ verb, noun, isOpen, ...p }) => (
  <Modal isOpen={isOpen} isCentered closeOnOverlayClick={false} size="xl">
    <ModalOverlay />
    <ModalContent>
      <ModalBody>
        <LoadingIndicator verb={verb} noun={noun} {...p} />
      </ModalBody>
    </ModalContent>
  </Modal>
);
