import { Heading, HStack, Text } from "@chakra-ui/react";
import {
  DraggableListField,
  FileUploadField,
  FormattedTextInput,
  FormikInput,
  ImageUploadField,
  NumberInputField,
  OptionsField,
  TextAreaField,
} from "components/core/forms";
// import {
//   ChemicalDisposableTable,
//   GreenMetricsCalculator,
//   GroupPlanTable,
//   HazardSummaryTable,
//   ProductYieldTable,
//   ReactionScheme,
//   TabbedImportPanel,
// } from "components/experiment-forms";
// import { Feedback } from "components/feedback/feedback";
import { FIELDS, INPUT_TYPES, SECTION_TYPES } from "constants";

const FieldWrapper = ({ hasFeedback = true, field, item, children }) => (
  <HStack>
    {children}
    {hasFeedback && <Feedback field={field} item={item} />}
  </HStack>
);

const Header = ({ field }) => (
  <Heading size="sm" as="u">
    {field.name}
  </Heading>
);

const Content = ({ field }) => (
  <HStack>
    <Text fontSize="sm" fontWeight="semibold">
      {field.name}
    </Text>
    <Text fontSize="sm">{field.defaultResponse}</Text>
  </HStack>
);

const TextField = ({ field, isDisabled }) => (
  <FormikInput
    name={field.id}
    label={field.name}
    isRequired={field.mandatory}
    placeholder={field.name}
    isDisabled={isDisabled}
  />
);

const FormattedTextField = ({ field, isDisabled }) => (
  <FormattedTextInput
    name={field.id}
    label={field.name}
    isRequired={field.mandatory}
    placeholder={field.name}
    isDisabled={isDisabled}
  />
);

const DateAndTimeField = ({ field, isDisabled }) => (
  <FormikInput
    name={field.id}
    label={field.name}
    isRequired={field.mandatory}
    placeholder={field.name}
    isDisabled={isDisabled}
    type="datetime-local"
  />
);

const NumberField = ({ field, isDisabled }) => (
  <NumberInputField
    name={field.id}
    label={field.name}
    isRequired={field.mandatory}
    placeholder={field.name}
    isDisabled={isDisabled}
  />
);

const DescriptionField = ({ field, isDisabled }) => (
  <TextAreaField
    name={field.id}
    title={field.name}
    placeholder={field.name}
    isRequired={field.mandatory}
    isDisabled={isDisabled}
  />
);

const FileField = ({ field, isDisabled }) => (
  <FileUploadField
    name={field.id}
    title={field.name}
    accept={field.fieldResponse?.accept ?? [".pdf", ".docx", ".doc"]}
    isRequired={field.mandatory}
    isDisabled={isDisabled}
  />
);

const ImageFileField = ({ field, isDisabled }) => (
  <ImageUploadField
    name={field.id}
    title={field.name}
    accept={field.fieldResponse?.accept ?? [".png", ".jpg", ".jpeg"]}
    isRequired={field.mandatory}
    isDisabled={isDisabled}
  />
);

const SortableListField = ({ field, isDisabled }) => (
  <DraggableListField
    name={field.id}
    label={field.name}
    isDisabled={isDisabled}
  />
);

const ReactionSchemeField = ({ field, isDisabled }) => (
  <ReactionScheme name={field.id} isDisabled={isDisabled} />
);

const MultiReactionSchemeField = ({ field, isDisabled }) => (
  <TabbedImportPanel
    name={field.id}
    label={field.name}
    isDisabled={isDisabled}
    sourceType={SECTION_TYPES.Note}
    fieldName={FIELDS.ReactionSchemeField}
    Component={ReactionScheme}
  />
);

const ChemicalDisposalTableField = ({ field, isDisabled }) => (
  <ChemicalDisposableTable
    name={field.id}
    label={field.name}
    isDisabled={isDisabled}
  />
);

const MultipleField = ({ field, isDisabled }) => (
  <OptionsField
    name={field.id}
    label={field.name}
    options={field.selectFieldOptions}
    isMultiple
    isDisabled={isDisabled}
    isRequired={field.mandatory}
  />
);

const RadioField = ({ field, isDisabled }) => (
  <OptionsField
    name={field.id}
    label={field.name}
    options={field.selectFieldOptions}
    isDisabled={isDisabled}
    isRequired={field.mandatory}
  />
);

const ProjectGroupPlanTableField = ({ field, isDisabled }) => (
  <GroupPlanTable name={field.id} label={field.name} isDisabled={isDisabled} />
);

const ProjectGroupHazardTableField = ({ field, isDisabled }) => (
  <HazardSummaryTable
    name={field.id}
    label={field.name}
    isDisabled={isDisabled}
  />
);

const YieldTableField = ({ field, isDisabled }) => (
  <ProductYieldTable
    name={field.id}
    label={field.name}
    isDisabled={isDisabled}
  />
);

const MultiYieldTableField = ({ field, isDisabled }) => (
  <TabbedImportPanel
    name={field.id}
    label={field.name}
    isDisabled={isDisabled}
    sourceType={SECTION_TYPES.Note}
    fieldName={FIELDS.YieldCalculationField}
    Component={ProductYieldTable}
  />
);

const GreenMetricsTableField = ({ field, isDisabled }) => (
  <GreenMetricsCalculator name={field.id} isDisabled={isDisabled} />
);

const MultiGreenMetricsTableField = ({ field, isDisabled }) => (
  <TabbedImportPanel
    name={field.id}
    label={field.name}
    isDisabled={isDisabled}
    sourceType={SECTION_TYPES.Note}
    fieldName={FIELDS.GreenMetricsTable}
    Component={GreenMetricsCalculator}
  />
);

export const INPUT_TYPES_MAP = {
  [INPUT_TYPES.Content]: Content,
  [INPUT_TYPES.Text]: TextField,
  [INPUT_TYPES.Description]: DescriptionField,
  [INPUT_TYPES.Number]: NumberField,
  [INPUT_TYPES.ReactionScheme]: ReactionSchemeField,
  [INPUT_TYPES.MultiReactionScheme]: MultiReactionSchemeField,
  [INPUT_TYPES.Radio]: RadioField,
  [INPUT_TYPES.Multiple]: MultipleField,
  [INPUT_TYPES.File]: FileField,
  [INPUT_TYPES.Header]: Header,
  [INPUT_TYPES.ChemicalDisposalTable]: ChemicalDisposalTableField,
  [INPUT_TYPES.SortableList]: SortableListField,
  [INPUT_TYPES.ProjectGroupPlanTable]: ProjectGroupPlanTableField,
  [INPUT_TYPES.ProjectGroupHazardTable]: ProjectGroupHazardTableField,
  [INPUT_TYPES.ImageFile]: ImageFileField,
  [INPUT_TYPES.YieldTable]: YieldTableField,
  [INPUT_TYPES.MultiYieldTable]: MultiYieldTableField,
  [INPUT_TYPES.GreenMetricsTable]: GreenMetricsTableField,
  [INPUT_TYPES.MultiGreenMetricsTable]: MultiGreenMetricsTableField,
  [INPUT_TYPES.DateAndTime]: DateAndTimeField,
  [INPUT_TYPES.FormattedTextInput]: FormattedTextField,
};

const NO_WRAPPER_COMPONENTS = [Header, Content];

/**
 * Creates a field based on its input type.
 */
export const Field = ({ field, isDisabled, item }) => {
  const inputType = field.inputType.name.toUpperCase();
  const [, Component] =
    Object.entries(INPUT_TYPES_MAP).find(
      ([key]) => key.toUpperCase() === inputType,
    ) || [];

  if (!Component) {
    return null;
  }

  return NO_WRAPPER_COMPONENTS.includes(Component) ? (
    <Component field={field} isDisabled={isDisabled} />
  ) : (
    <FieldWrapper field={field} item={item}>
      <Component field={field} isDisabled={isDisabled} />
    </FieldWrapper>
  );
};
