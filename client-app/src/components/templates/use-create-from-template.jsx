import { useBackendApi } from "contexts";
import { useToast } from "@chakra-ui/react";
import { useRef, useState } from "react";

export const useCreateFormFromTemplate = ({
  templateType,
  reactionId,
  workgroupName,
  workbookName,
  onCreated,
}) => {
  const { projects: api } = useBackendApi();
  const toast = useToast();
  const formRef = useRef();

  const [isLoading, setIsLoading] = useState(false);
  const [feedback, setFeedback] = useState(null);

  const initialValues = { templateId: [] };

  const handleSubmit = async (values) => {
    try {
      setIsLoading(true);
      setFeedback(null);

      const response = await api.create({
        reactionId,
        workgroupName,
        workbookName,
        templateId: Number(values.templateId[0]),
        templateType: templateType,
      });

      if (!response.ok) throw new Error("Failed to create form");

      const data = await response.json();

      if (data?.id || data?.uuid) {
        toast({
          title: templateType + " created",
          status: "success",
          duration: 10000,
          isClosable: true,
          position: "top",
        });
        onCreated(data);
      }
    } catch (e) {
      console.error(e);
      setFeedback({
        status: "error",
        message: "Failed to create " + templateType,
      });
    } finally {
      setIsLoading(false);
    }
  };

  return {
    formRef,
    initialValues,
    feedback,
    isLoading,
    handleSubmit,
  };
};
