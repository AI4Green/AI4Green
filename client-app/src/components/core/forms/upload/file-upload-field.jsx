import {
  Alert,
  AlertDescription,
  AlertIcon,
  Button,
  FormControl,
  FormLabel,
  HStack,
  Input,
  Link,
  Tag,
  TagCloseButton,
  TagLabel,
  TagLeftIcon,
  Text,
  useDisclosure,
  VStack,
} from "@chakra-ui/react";
import { FormHelpError } from "components/core/forms";
import { Modal } from "components/core/modal.tsx";
import { useField } from "formik";
import { useRef, useState } from "react";
import { FaCloudUploadAlt } from "react-icons/fa";
import { MdCheckCircle } from "react-icons/md";

/**
 * This is a formik field component for uploading files.
 *
 * The structure of the field value is an array of objects:
 * - file: the file object
 * - name: the file name
 * - isNew: a flag to indicate if the file is newly uploaded
 * - isMarkedForDeletion: a flag to indicate if the file is marked for deletion (only applied to existing files)
 * - additional metadata can be added
 *
 */

const InfoAlert = ({ accept }) => {
  return (
    <Alert borderRadius={7} variant="left-accent" colorScheme="gray" py={2}>
      <AlertIcon />
      <HStack align="start">
        <Text fontSize="xs">Supported format</Text>
        <HStack>
          {accept?.map((extension, index) => (
            <Tag key={index} colorScheme="green">
              <TagLeftIcon boxSize="12px" as={MdCheckCircle} />
              <TagLabel>{extension}</TagLabel>
            </Tag>
          ))}
        </HStack>
      </HStack>
    </Alert>
  );
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

export const FileUploadField = ({
  name,
  title,
  accept,
  isRequired,
  isDisabled,
  isMultiple, // allow multiple files to be uploaded
}) => {
  const [field, meta, helpers] = useField(name);

  const [fileUploadError, setFileUploadError] = useState();

  const fileInputRef = useRef();

  const handleBtnClick = () => fileInputRef.current.click();

  const handleChange = ({ target: { files } }) => {
    const existingFiles = field.value;
    const selectedFile = files[0];
    const selectedFileExtension = selectedFile.name
      .slice(selectedFile.name.lastIndexOf("."))
      .toLowerCase();

    const isFileExtensionAccepted = accept.some(
      (extension) => extension.toLowerCase() === selectedFileExtension,
    );
    const isDuplicateName = existingFiles.some(
      (existingFile) => existingFile.name === selectedFile.name,
    );

    const fileData = {
      file: {
        file: isFileExtensionAccepted && !isDuplicateName ? selectedFile : null,
        name: selectedFile.name,
        isNew: true,
        // additional metadata goes here
      },
      error: !isFileExtensionAccepted
        ? `File ${selectedFile.name} is not supported!`
        : isDuplicateName
        ? `File ${selectedFile.name} already exists!`
        : null,
    };

    setFileUploadError(fileData.error);
    if (fileData.error) return;

    helpers.setValue([...existingFiles, fileData.file]);
  };

  const handleRemove = (f) => {
    const files = field.value.filter((file) => file.name !== f.name);
    helpers.setValue(files);
  };

  const handleRemoveExisting = (f) => {
    const updatedFiles = field.value.map((file) =>
      file.name === f.name ? { ...file, isMarkedForDeletion: true } : file,
    );
    helpers.setValue(updatedFiles);
  };

  const handleDownload = async (download, fileName) => {
    const blob = await download();
    const url = URL.createObjectURL(blob);
    const link = Object.assign(document.createElement("a"), {
      href: url,
      download: fileName,
    });
    document.body.appendChild(link).click();
    URL.revokeObjectURL(url);
  };

  const canUpload =
    !field.value?.length || // if no files uploaded
    (field.value.length === 1 && field.value[0].isMarkedForDeletion) || // and it's marked for deletion
    isMultiple; // or multiple files are allowed

  return (
    <FormControl id={field.name} isRequired={isRequired}>
      <VStack align="start" w="100%" spacing={2}>
        <FormLabel>{title}</FormLabel>

        {!isDisabled && (
          <>
            {canUpload && <InfoAlert accept={accept} />}
            {fileUploadError && (
              <FileUploadAlert status="error" message={fileUploadError} />
            )}

            {field.value
              ?.filter((file) => file.isNew)
              .map((f) => (
                <FileUploadAlert
                  key={f.name}
                  status="success"
                  message={`File ${f.name} successfully added!`}
                />
              ))}

            <HStack w="100%">
              {canUpload && ( // or multiple files are allowed
                <>
                  <Button
                    colorScheme="gray"
                    variant="outline"
                    size="sm"
                    leftIcon={<FaCloudUploadAlt />}
                    onClick={handleBtnClick}
                  >
                    {isMultiple ? "Select files" : "Select file"}
                  </Button>
                  <Input
                    name={name}
                    type="file"
                    display="none"
                    accept={accept}
                    onChange={handleChange}
                    ref={fileInputRef}
                  />
                </>
              )}

              <VStack align="start" w="100%">
                {field.value
                  ?.filter((file) => file.isNew)
                  .map((f) => (
                    <FileName
                      key={f.name}
                      fileName={f.name}
                      handleRemove={() => handleRemove(f)}
                    />
                  ))}
              </VStack>
            </HStack>
          </>
        )}

        {field.value
          ?.filter((file) => !file.isMarkedForDeletion && !file.isNew)
          .map((f) => (
            <Tag key={f.name} colorScheme="green">
              <TagLabel>
                <Link
                  onClick={async () => await handleDownload(f.download, f.name)}
                >
                  {f.name}
                </Link>
              </TagLabel>

              {!isDisabled && (
                <RemoveExisitngButton
                  file={f}
                  handleRemoveExisting={handleRemoveExisting}
                />
              )}
            </Tag>
          ))}
      </VStack>
      <FormHelpError
        isInvalid={meta.touched && meta.error}
        error={meta.error}
        collapseEmpty
        replaceHelpWithError
      />
    </FormControl>
  );
};

const RemoveExisitngButton = ({ file, handleRemoveExisting }) => {
  const { isOpen, onOpen, onClose } = useDisclosure();
  return (
    <>
      <TagCloseButton onClick={onOpen} />
      {isOpen && (
        <Modal
          body={modalBody}
          title="🚫 File removal warning"
          actionBtnCaption="Continue"
          actionBtnColorScheme="orange"
          onAction={() => handleRemoveExisting(file)}
          isOpen={isOpen}
          onClose={onClose}
        />
      )}
    </>
  );
};

const FileName = ({ fileName, handleRemove }) => (
  <Tag variant="outline" py={1}>
    <TagLabel>{fileName}</TagLabel>
    <TagCloseButton onClick={handleRemove} color="red" />
  </Tag>
);

const FileUploadAlert = ({ status, message }) => {
  return (
    <Alert status={status} py={2}>
      <AlertIcon />
      <AlertDescription>
        <Text fontSize="xs">{message}</Text>
      </AlertDescription>
    </Alert>
  );
};
