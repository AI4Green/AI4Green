import {
  Box,
  Button,
  HStack,
  Image,
  Input,
  Tag,
  TagCloseButton,
  TagLabel,
  TagLeftIcon,
  VStack,
} from "@chakra-ui/react";
import { useObjectUrl } from "helpers/hooks";
import { useRef } from "react";
import { FaCheckCircle, FaCloudUploadAlt } from "react-icons/fa";

import { CaptionInput } from "./image-upload-field";

/**
 * Displays an upload button that opens the dialog to select images.
 */
export const UploadImage = ({
  accept,
  setFileUploadError,
  values,
  helpers,
}) => {
  const fileInputRef = useRef();
  const handleClick = () => fileInputRef.current.click();

  const handleChange = ({ target: { files } }) => {
    const existingImages = values;
    let errors = [];

    const newFiles = Array.from(files).map((selectedFile) => {
      const selectedFileExtension = selectedFile.name
        .slice(selectedFile.name.lastIndexOf("."))
        .toLowerCase();

      const isFileExtensionAccepted = accept.some(
        (extension) => extension.toLowerCase() === selectedFileExtension,
      );
      const isDuplicateName = existingImages.some(
        (existingFile) => existingFile.name === selectedFile.name,
      );

      if (!isFileExtensionAccepted) {
        errors.push({
          fileName: selectedFile.name,
          message: `File ${selectedFile.name} is not supported!`,
        });
        return null;
      }

      if (isDuplicateName) {
        errors.push({
          fileName: selectedFile.name,
          message: `File ${selectedFile.name} already exists!`,
        });
        return null;
      }

      return {
        file: selectedFile,
        name: selectedFile.name,
        isNew: true,
        caption: selectedFile.name,
      };
    });

    setFileUploadError(errors);
    helpers.setValue([...existingImages, ...newFiles.filter(Boolean)]);
  };
  return (
    <Box>
      <Button
        colorScheme="gray"
        variant="outline"
        size="sm"
        leftIcon={<FaCloudUploadAlt />}
        onClick={handleClick}
      >
        Select Images
      </Button>
      <Input
        type="file"
        display="none"
        accept={accept}
        onChange={handleChange}
        ref={fileInputRef}
        multiple
      />
    </Box>
  );
};

/**
 * Displays an image that has been selected for upload with caption input field and remove button.
 */
export const ImageForUpload = ({
  setFileUploadError,
  fileUploadError,
  image,
  values,
  helpers,
}) => {
  const objectUrl = useObjectUrl(image.file);

  // Remove the image from the list of images to be uploaded
  const handleRemove = () => {
    const files = values.filter((item) => item.name !== image.name);
    setFileUploadError(
      fileUploadError.filter((error) => error.fileName !== image.name),
    );
    helpers.setValue(files);
  };

  return (
    <Box>
      <HStack>
        <Image
          src={objectUrl}
          alt={image.name}
          h={24}
          maxW={28}
          objectFit="contain"
        />
        <VStack align="start" spacing={2}>
          <CaptionInput image={image} values={values} helpers={helpers} />

          <Tag variant="outline" size="sm" borderRadius="full">
            <TagLeftIcon as={FaCheckCircle} color="green" />
            <TagLabel>{image.name}</TagLabel>
            <TagCloseButton onClick={handleRemove} color="red" />
          </Tag>
        </VStack>
      </HStack>
    </Box>
  );
};
