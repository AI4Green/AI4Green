import { useCoshhFormsList } from "api";
import { useCreateFormFromTemplate } from "components/templates";
import { CreateFromTemplateInline } from "components/templates/create.jsx";

export const CoshhCreateInline = ({
  reactionId,
  workgroupName,
  workbookName,
  onCreated,
  isOpen,
  onClose,
}) => {
  const { data: templates = [] } = useCoshhFormsList();

  const { formRef, initialValues, feedback, isLoading, handleSubmit } =
    useCreateFormFromTemplate({
      templateType: "COSHH",
      reactionId,
      workgroupName,
      workbookName,
      onCreated,
    });

  return (
    <CreateFromTemplateInline
      title="Create COSHH Form"
      templateLabel="COSHH"
      formRef={formRef}
      templates={templates}
      initialValues={initialValues}
      feedback={feedback}
      isLoading={isLoading}
      onSubmit={handleSubmit}
      isOpen={isOpen}
      onClose={onClose}
    />
  );
};
