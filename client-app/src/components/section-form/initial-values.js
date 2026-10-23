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
