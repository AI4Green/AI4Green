import { INPUT_TYPES } from "constants";

const parseFileValue = (raw) => {
  if (!raw) return [];
  try {
    const parsed = typeof raw === "string" ? JSON.parse(raw) : raw;
    return Array.isArray(parsed) ? parsed : [parsed];
  } catch {
    return [];
  }
};

const getInitialValue = (field, recordId, sectionId, sectionApis) => {
  const {
    File,
    ImageFile,
    Text,
    Description,
    Multiple,
    Radio,
    SortableList,
    ChemicalDisposalTable,
    Header,
    Content,
    ProjectGroupPlanTable,
    ProjectGroupHazardTable,
    MultiReactionScheme,
    MultiYieldTable,
    MultiGreenMetricsTable,
  } = INPUT_TYPES;

  const fieldType = field.inputType.title.toUpperCase();
  const fieldResponse = field.response?.value ?? "";

  switch (fieldType) {
    case File.toUpperCase():
    case ImageFile.toUpperCase(): {
      return {
        [field.id]: parseFileValue(fieldResponse).map((meta) => ({
          name: meta.display_name,
          uuid: meta.uuid,
          instanceId: recordId, // plain data, so the download call knows which instance
          isNew: false,
        })),
      };
    }

    case Text.toUpperCase():
    case Description.toUpperCase():
      return { [field.id]: fieldResponse ?? "" };

    case Multiple.toUpperCase():
    case Radio.toUpperCase():
    case SortableList.toUpperCase():
    case ChemicalDisposalTable.toUpperCase():
    case ProjectGroupPlanTable.toUpperCase():
    case ProjectGroupHazardTable.toUpperCase():
    case MultiReactionScheme.toUpperCase():
    case MultiYieldTable.toUpperCase():
    case MultiGreenMetricsTable.toUpperCase():
      return {
        [field.id]: !fieldResponse ? [] : fieldResponse,
      };

    case Header.toUpperCase():
    case Content.toUpperCase():
      return {};

    default:
      return { [field.id]: fieldResponse };
  }
};

/**
 * Gets a field's saved response: the response record's id and its most recent
 * value. Uses `??` rather than `||` so falsy answers such as `0` or `false` are
 * kept instead of being treated as unanswered.
 * @param {object} field - A template field, including its `fieldResponse` list.
 * @returns {{ id: (number|null), value: * }} The response id and latest value,
 *   each `null` when there is no response.
 */
export const getLatestResponse = (field) => {
  const response = field.fieldResponse?.[0];
  return {
    id: response?.id ?? null,
    value: response?.fieldResponseValues?.at(-1)?.value ?? null,
  };
};

export const initialValues = (
  sectionFields,
  recordId,
  sectionId,
  sectionApis,
) =>
  sectionFields.reduce((acc, field) => {
    return {
      ...acc,
      ...getInitialValue(field, recordId, sectionId, sectionApis),
    };
  }, {});
