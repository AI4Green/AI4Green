import {
  Box,
  Flex,
  HStack,
  IconButton,
  Image,
  Text,
  useDisclosure,
  VStack,
} from "@chakra-ui/react";
import { Modal } from "components/core/modal";
import { useObjectUrl } from "helpers/hooks";
import { useEffect, useState } from "react";
import { FaDownload, FaTimes } from "react-icons/fa";

import { CaptionInput } from "./image-upload-field";

/**
 * Displays an uploaded image with caption and download and remove (mark image for deletion) buttons.
 */
export const UploadedImage = ({ values, helpers, image }) => {
  const [blob, setBlob] = useState(null);
  const objectUrl = useObjectUrl(blob);

  useEffect(() => {
    image.download().then(setBlob);
  }, [image]);

  const handleDownload = async () => {
    const { name, caption, download } = image;
    const blob = await download();
    const url = URL.createObjectURL(blob);
    const link = Object.assign(document.createElement("a"), {
      href: url,
      download: `${caption}_${name}`,
    });
    document.body.appendChild(link).click();
    URL.revokeObjectURL(url);
  };

  return (
    <Flex
      w={48}
      borderWidth={1}
      borderRadius={4}
      overflow="hidden"
      direction="column"
      justify="space-between"
    >
      <Box>
        <Image src={objectUrl} h={24} maxW="xs" objectFit="contain" />
        <CaptionInput image={image} values={values} helpers={helpers} />
      </Box>
      <HStack justify="flex-end" mt="4">
        <RemoveExisitngButton image={image} values={values} helpers={helpers} />
        <IconButton
          aria-label="Download file"
          icon={<FaDownload />}
          size="xs"
          variant="ghost"
          onClick={handleDownload}
          colorScheme="green"
        />
      </HStack>
    </Flex>
  );
};

// Button to remove (mark for deletion) an image that was uploaded previously
const RemoveExisitngButton = ({ image, values, helpers }) => {
  const { isOpen, onOpen, onClose } = useDisclosure();

  // Mark the image for deletion (only applies to images that were uploaded previously)
  // Once the form is saved, the image will be deleted from the server
  const handleRemoveExisting = () => {
    const updatedFiles = values.map((item) =>
      item.name === image.name ? { ...item, isMarkedForDeletion: true } : item,
    );
    helpers.setValue(updatedFiles);
  };

  const modalBody = (
    <VStack align="flex-start">
      <Text fontSize="md" color="gray.600">
        You have chosen to remove the uploaded file.
      </Text>
      <Text fontSize="sm">
        By proceeding and saving the form, your previously uploaded file will be
        permanently deleted and cannot be recovered.
      </Text>
    </VStack>
  );
  return (
    <>
      <IconButton
        aria-label="Remove image"
        icon={<FaTimes />}
        size="xs"
        variant="ghost"
        onClick={onOpen}
        colorScheme="red"
      />
      {isOpen && (
        <Modal
          body={modalBody}
          title="🚫 File removal warning"
          actionBtnCaption="Continue"
          actionBtnColorScheme="orange"
          onAction={handleRemoveExisting}
          isOpen={isOpen}
          onClose={onClose}
        />
      )}
    </>
  );
};
