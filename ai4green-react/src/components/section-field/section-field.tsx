import { SECTION_TYPES, STAGES_PERMISSIONS } from "constants";
import { useMemo } from "react";

import { Field, isFieldTriggered } from ".";

/**
 * Creates a fields for the section form.
 * Mainly used to create a input field (formik field).
 * Also, renders a feedback component next to the input field.
 * Feedback component allows the instructor to approve or leave a comment on the field response.
 * Similarly, it allows the student view the feedback.
 * @param {*}
 * field:
 *  - field object mainly used for creating the field.
 *  - contains properties like id, name, fieldType, mandatory, trigger, etc. used to create the field.
 * fieldValues: collection of field values, which can be accessed by using the field.id as key
 * recordId: could be planId or reportId
 * sectionFields: collection of fields in the section
 * @returns
 */
export const SectionField = ({
  item,
  field,
  fields,
  fieldValues,
  feedback,
  isInstructor,
}) => {
  const { canEdit } = useFieldPermissions(item, field, isInstructor);

  return (
    <>
      <Field field={{ ...field, feedback }} isDisabled={!canEdit} item={item} />
      {field.triggerField && (
        <TriggerField
          item={item}
          field={field}
          fields={fields}
          fieldValues={fieldValues}
          isInstructor={isInstructor}
        />
      )}
    </>
  );
};

const TriggerField = ({ item, field, fields, fieldValues, isInstructor }) => {
  const { isFieldTriggeringChild, triggerTargetField } = useTriggerField(
    field,
    fields,
    fieldValues,
  );

  if (!isFieldTriggeringChild || !triggerTargetField) {
    return null;
  }

  return (
    <>
      <SectionField
        item={item}
        field={triggerTargetField}
        fields={fields}
        fieldValues={fieldValues}
        feedback={triggerTargetField.feedback}
        isInstructor={isInstructor}
      />
      {triggerTargetField?.triggerField && (
        <TriggerField
          item={item}
          field={triggerTargetField}
          fields={fields}
          fieldValues={fieldValues}
          isInstructor={isInstructor}
        />
      )}
    </>
  );
};

const useFieldPermissions = (item, field, isInstructor) => {
  const hasRequiredPermissions = [
    STAGES_PERMISSIONS.OwnerCanEdit,
    STAGES_PERMISSIONS.OwnerCanEditCommented,
  ].some(
    (permission) =>
      item.stage?.permissions.includes(permission) && !field.feedback?.approved,
  );

  const canEdit =
    !isInstructor &&
    (item.type.toUpperCase() === SECTION_TYPES.ProjectGroup.toUpperCase() ||
      (item.isOwner && hasRequiredPermissions));

  return { canEdit };
};

const useTriggerField = (field, fields, fieldValues) => {
  const isFieldTriggeringChild = isFieldTriggered(
    field.inputType.name,
    field.triggerField?.value,
    fieldValues[field.id],
  );

  // get the trigger target field from the fields collection
  const triggerTargetField = useMemo(
    () => fields.find((x) => x.id === field.triggerField?.id),
    [fields, field.triggerField?.id],
  );

  return { isFieldTriggeringChild, triggerTargetField };
};
