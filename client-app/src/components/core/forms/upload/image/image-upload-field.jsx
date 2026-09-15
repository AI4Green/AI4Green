import {
  Alert,
  AlertDescription,
  AlertIcon,
  Flex,
  FormControl,
  FormLabel,
  HStack,
  Input,
  Tag,
  TagLabel,
  TagLeftIcon,
  Text,
  VStack,
} from "@chakra-ui/react";
import { FormHelpError } from "components/core/forms";
import { useField } from "formik";
import { useMemo, useState } from "react";
import { MdCheckCircle } from "react-icons/md";

import { ImageForUpload, UploadImage } from "./upload";
import { UploadedImage } from "./uploaded";

/**
 * This is a formik field component for uploading images with captions.
 *
 * The structure of the field value is an array of objects:
 * - file: the file object
 * - name: the file name
 * - isNew: a flag to indicate if the file is newly uploaded
 * - isMarkedForDeletion: a flag to indicate if the file is marked for deletion (only applied to existing files)
 * - caption: the caption of the image
 *
 */

export const ImageUploadField = ({
  name,
  title,
  accept,
  isRequired,
  isDisabled,
}) => {
  const [field, meta, helpers] = useField(name);
  const [fileUploadError, setFileUploadError] = useState([]);

  const processedFieldValue = useMemo(() => {
    return field.value?.map((file, i) => {
      if (Array.isArray(meta.error)) {
        return { ...file, error: { caption: meta.error[i]?.caption } };
      }
      const { error, ...other } = file; // remove error key when there is no error
      return other;
    });
  }, [field.value, meta.error]);

  return (
    <FormControl
      id={field.name}
      isRequired={isRequired}
      isInvalid={meta.error && meta.touched}
    >
      <VStack align="start" w="100%" spacing={2}>
        <FormLabel>{title}</FormLabel>

        {!isDisabled && (
          <>
            <Info accept={accept} />
            {fileUploadError &&
              fileUploadError.map((e, i) => (
                <ImageForUploadAlert
                  key={i}
                  status="error"
                  message={e.message}
                />
              ))}

            <VStack w="full" align="flex-start">
              <UploadImage
                accept={accept}
                setFileUploadError={setFileUploadError}
                values={processedFieldValue}
                helpers={helpers}
              />

              <VStack align="start" w="100%">
                {processedFieldValue
                  ?.filter((image) => image.isNew)
                  .map((item) => (
                    <ImageForUpload
                      key={item.name}
                      setFileUploadError={setFileUploadError}
                      fileUploadError={fileUploadError}
                      image={item}
                      values={processedFieldValue}
                      helpers={helpers}
                    />
                  ))}
              </VStack>
            </VStack>
          </>
        )}
        <Flex wrap="wrap" gap={2}>
          {processedFieldValue
            ?.filter((file) => !file.isMarkedForDeletion && !file.isNew)
            .map((item) => (
              <UploadedImage
                key={item.name}
                image={item}
                values={processedFieldValue}
                helpers={helpers}
              />
            ))}
        </Flex>
      </VStack>
      <FormHelpError
        isInvalid={meta.touched && meta.error && !Array.isArray(meta.error)} // only show if not an array
        error={meta.error}
        collapseEmpty
        replaceHelpWithError
      />
    </FormControl>
  );
};

const Info = ({ accept }) => {
  return (
    <Alert borderRadius={7} variant="left-accent" colorScheme="gray" py={2}>
      <AlertIcon />
      <HStack align="start">
        <Text fontSize="xs">Supported format</Text>
        <HStack>
          {accept?.map((extension, index) => (
            <Tag key={index} variant="subtle" colorScheme="green">
              <TagLeftIcon boxSize="12px" as={MdCheckCircle} />
              <TagLabel>{extension}</TagLabel>
            </Tag>
          ))}
        </HStack>
      </HStack>
    </Alert>
  );
};

// Displays whether or not file has been added to the list of files to be uploaded
const ImageForUploadAlert = ({ status, message }) => {
  return (
    <Alert status={status} py={2}>
      <AlertIcon boxSize={4} />
      <AlertDescription>
        <Text fontSize="xs">{message}</Text>
      </AlertDescription>
    </Alert>
  );
};

/**
 * Input field for inputting image caption.
 */
export const CaptionInput = ({ image, values, helpers }) => {
  // Handle the image caption input changes
  const handleCaptionInputChange = (e) => {
    const updatedImages = [...values];
    const fileIndex = updatedImages.findIndex(
      (item) => item.name === image.name,
    );
    updatedImages[fileIndex] = { ...image, caption: e.target.value };
    helpers.setValue(updatedImages);
  };
  return (
    <FormControl isInvalid={image.error?.caption}>
      <Input
        mt={4}
        px={2}
        variant="flushed"
        placeholder="Add image caption"
        size="xs"
        onChange={handleCaptionInputChange}
        value={image?.caption}
        color="gray.500"
      />
      <FormHelpError
        isInvalid={image.error?.caption}
        error={image.error?.caption}
        collapseEmpty
        replaceHelpWithError
      />
    </FormControl>
  );
};
