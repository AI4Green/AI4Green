import { INPUT_TYPES } from "constants";

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

  const fieldType = field.inputType.name.toUpperCase();
  const fieldResponse = field.response.value;

  switch (fieldType) {
    case File.toUpperCase():
    case ImageFile.toUpperCase(): {
      return {
        [field.id]: Array.isArray(fieldResponse)
          ? fieldResponse.map((file) => ({
              ...file,
              download: async () => {
                const response = await sectionApis.downloadSectionFile(
                  sectionId,
                  recordId,
                  file.location,
                  file.name,
                );
                return await response.blob();
              },
            }))
          : [],
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
