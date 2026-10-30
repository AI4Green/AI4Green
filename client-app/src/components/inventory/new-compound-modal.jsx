import { Grid, Text, VStack } from "@chakra-ui/react";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { number, object, string } from "yup";
import { SelectField, TextField } from "./form-fields.jsx";
import { ErrorMessage, InvModal, MOBILE } from "./ui.jsx";
import { useInventoryAction } from "./use-inventory-action.js";
import { CAS_PATTERN } from "./utils.jsx";

const FORM_ID = "inv-new-compound-form";

// dash-separated GHS codes, e.g. H225-H319 or EUH066, as the reaction table uses
const HAZARD_CODES_PATTERN =
  /^(EU)?H\d{3}[A-Za-z]{0,2}(-(EU)?H\d{3}[A-Za-z]{0,2})*$/;

/**
 * Optional positive number field: blank is allowed, anything else must be a
 * number between 0 and `max`. The upper bound also rejects Infinity.
 * @param {number} max The largest value allowed.
 */
const optionalPositive = (max) =>
  number()
    .transform((value, original) => (original === "" ? undefined : value))
    .typeError("Enter a number")
    .positive("Must be more than zero")
    .max(max, `Must be ${max} or less`);

const schema = object({
  name: string()
    .trim()
    .required("Enter the compound's name")
    .max(200, "Keep the name under 200 characters"),
  cas: string()
    .trim()
    .required("Enter a CAS number")
    .matches(CAS_PATTERN, "Enter a CAS number like 67-64-1"),
  workbookId: string().required("Choose a workbook"),
  smiles: string().trim().max(1000, "Keep the SMILES under 1000 characters"),
  molWeight: optionalPositive(100000),
  density: optionalPositive(100),
  concentration: optionalPositive(1000),
  hazardCodes: string().trim().matches(HAZARD_CODES_PATTERN, {
    message: "Use H-codes separated by dashes, e.g. H225-H319",
    excludeEmptyString: true,
  }),
});

/** Turns an optional numeric field into a number, or null when blank. */
const optionalNumber = (value) => (value === "" ? null : Number(value));

/**
 * Add a compound that isn't in the compound database, as a novel compound in
 * one of the workgroup's workbooks (the same place the reaction table puts
 * them, so it can be used in that workbook's reactions too). Opened from the
 * compound lookup when a CAS number has no match; shown on top of Add stock.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup adding the compound.
 * @param {object[]} props.workbooks The workgroup's workbooks, as `{ id, name }`.
 * @param {string} props.cas The CAS number that wasn't found, to start from.
 * @param {(compound: object) => void} props.onCreated Called with the new compound.
 * @param {() => void} props.onClose Closes the modal.
 */
export const NewCompoundModal = ({
  workgroupName,
  workbooks,
  cas,
  onCreated,
  onClose,
}) => {
  const { inventory } = useBackendApi();
  const { run, isLoading, error } = useInventoryAction(workgroupName);

  const handleSubmit = async (values) => {
    let created = null;
    const ok = await run(async () => {
      created = await inventory.addNovelCompound(workgroupName, {
        workbookId: Number(values.workbookId),
        name: values.name.trim(),
        cas: values.cas.trim(),
        smiles: values.smiles.trim() || null,
        molWeight: optionalNumber(values.molWeight),
        density: optionalNumber(values.density),
        concentration: optionalNumber(values.concentration),
        hazardCodes: values.hazardCodes.trim() || null,
      });
    }, "Compound added");
    if (ok) {
      onCreated(created);
      onClose();
    }
  };

  return (
    <InvModal
      isOpen
      onClose={onClose}
      title="Add new compound"
      actionLabel="Add compound"
      formId={FORM_ID}
      isLoading={isLoading}
      size="xl"
    >
      {error && <ErrorMessage>{error}</ErrorMessage>}
      <Text mb="16px">
        This CAS number isn't in the compound database. Add it as a novel
        compound in one of your workbooks, the same way as in the reaction
        table.
      </Text>
      <Formik
        initialValues={{
          name: "",
          cas,
          workbookId: workbooks.length === 1 ? String(workbooks[0].id) : "",
          smiles: "",
          molWeight: "",
          density: "",
          concentration: "",
          hazardCodes: "",
        }}
        validationSchema={schema}
        onSubmit={handleSubmit}
      >
        <VStack as={Form} id={FORM_ID} noValidate spacing="16px">
          <TextField name="name" label="Name" isRequired autoFocus />
          <Grid
            templateColumns="1fr 1fr"
            gap="16px"
            w="100%"
            sx={{ [MOBILE]: { gridTemplateColumns: "1fr" } }}
          >
            <TextField name="cas" label="CAS number" isRequired />
            <SelectField
              name="workbookId"
              label="Workbook"
              isRequired
              helpText={
                workbooks.length === 0
                  ? "This workgroup has no workbooks yet. Create one first."
                  : "Novel compounds belong to a workbook."
              }
            >
              <option value="">Choose a workbook</option>
              {workbooks.map((wb) => (
                <option key={wb.id} value={wb.id}>
                  {wb.name}
                </option>
              ))}
            </SelectField>
          </Grid>
          <TextField
            name="smiles"
            label="SMILES"
            helpText="Used to work out the molecular weight if you leave it blank."
          />
          <Grid
            templateColumns="repeat(3, 1fr)"
            gap="16px"
            w="100%"
            sx={{ [MOBILE]: { gridTemplateColumns: "1fr" } }}
          >
            <TextField
              name="molWeight"
              label="Molecular weight (g/mol)"
              type="number"
              inputMode="decimal"
              step="any"
              min={0}
            />
            <TextField
              name="density"
              label="Density (g/mL)"
              type="number"
              inputMode="decimal"
              step="any"
              min={0}
            />
            <TextField
              name="concentration"
              label="Concentration (mol/L)"
              type="number"
              inputMode="decimal"
              step="any"
              min={0}
            />
          </Grid>
          <TextField
            name="hazardCodes"
            label="Hazard codes"
            placeholder="e.g. H225-H319"
            helpText="Separate codes with dashes. Leave blank if unknown."
          />
        </VStack>
      </Formik>
    </InvModal>
  );
};
