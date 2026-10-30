import { Grid, VStack } from "@chakra-ui/react";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { useRef, useState } from "react";
import { mixed, number, object, string } from "yup";
import { CompoundLookup, heldCompound } from "./compound-lookup.jsx";
import { SelectField, TextField } from "./form-fields.jsx";
import { NewCompoundModal } from "./new-compound-modal.jsx";
import { ErrorMessage, InvModal, MOBILE } from "./ui.jsx";
import { useInventoryAction } from "./use-inventory-action.js";
import { todayIso } from "./utils.jsx";

const UNITS = ["mg", "g", "kg", "µL", "mL", "L"];

const FORM_ID = "inv-add-stock-form";

/** Two fields side by side, stacked on narrow screens. */
const FieldPair = (props) => (
  <Grid
    templateColumns="1fr 1fr"
    gap="16px"
    w="100%"
    sx={{ [MOBILE]: { gridTemplateColumns: "1fr" } }}
    {...props}
  />
);

const schema = object({
  // the chosen compound object, set by CompoundLookup
  compound: mixed().test(
    "chosen",
    "Search for a compound and choose one",
    (value) => Boolean(value),
  ),
  locationId: string().required("Choose a location"),
  quantity: number()
    .typeError("Enter a number")
    .required("Enter the quantity")
    .moreThan(0, "Quantity must be more than zero"),
  unit: string().required("Choose a unit"),
  supplier: string().trim().max(255),
  lotNumber: string().trim().max(100),
  purity: string().trim().max(50),
  barcode: string().trim().max(100),
  dateReceived: string(),
  expiryDate: string(),
});

/**
 * Add a container of stock. With `workgroupCompound` set ("Add container" on
 * a compound row) the compound is fixed; without it ("+ Add stock") the user
 * looks one up by name or CAS number, and can add a compound that isn't in
 * the database as a novel compound.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup to add stock to.
 * @param {object[]} props.compounds The workgroup's WorkgroupCompounds.
 * @param {object[]} props.locations The workgroup's locations.
 * @param {object[]} props.workbooks The workgroup's workbooks, as `{ id, name }`.
 * @param {object} [props.workgroupCompound] The compound to add a container of.
 * @param {() => void} props.onClose Closes the modal.
 */
export const AddStockModal = ({
  workgroupName,
  compounds,
  locations,
  workbooks,
  workgroupCompound,
  onClose,
}) => {
  const { inventory } = useBackendApi();
  const { run, isLoading, error } = useInventoryAction(workgroupName);
  const formRef = useRef();
  // CAS number being added as a new compound, while that modal is open
  const [newCompoundCas, setNewCompoundCas] = useState(null);

  const activeLocations = locations.filter((l) => !l.isArchived);
  const isFixed = Boolean(workgroupCompound);

  const initialValues = {
    compound: isFixed ? heldCompound(workgroupCompound) : null,
    locationId: "",
    quantity: "",
    unit: workgroupCompound?.thresholdUnit ?? "",
    supplier: workgroupCompound?.preferredSupplier ?? "",
    lotNumber: "",
    purity: "",
    barcode: "",
    dateReceived: todayIso(),
    expiryDate: "",
  };

  const handleSubmit = async (values) => {
    const { compound } = values;
    const isHeld = Boolean(compound.workgroupCompoundId);
    const optional = (value) => value.trim() || null;

    const ok = await run(
      () =>
        inventory.addItem(workgroupName, {
          workgroupCompoundId: isHeld ? compound.workgroupCompoundId : null,
          // a compound the workgroup doesn't hold yet, from the lookup
          compound: isHeld
            ? null
            : { source: compound.source, id: compound.id },
          locationId: Number(values.locationId),
          quantity: Number(values.quantity),
          unit: values.unit,
          supplier: optional(values.supplier),
          lotNumber: optional(values.lotNumber),
          purity: optional(values.purity),
          barcode: optional(values.barcode),
          dateReceived: values.dateReceived || null,
          expiryDate: values.expiryDate || null,
        }),
      "Stock added",
    );
    if (ok) onClose();
  };

  // select the compound just created in the new compound modal
  const handleCompoundCreated = (compound) =>
    formRef.current?.setFieldValue("compound", compound);

  return (
    <InvModal
      isOpen
      onClose={onClose}
      title={
        isFixed
          ? `Add container of ${workgroupCompound.compound.name}`
          : "Add stock"
      }
      actionLabel="Add"
      formId={FORM_ID}
      isLoading={isLoading}
      size="xl"
    >
      {error && <ErrorMessage>{error}</ErrorMessage>}
      <Formik
        innerRef={formRef}
        initialValues={initialValues}
        validationSchema={schema}
        onSubmit={handleSubmit}
      >
        <VStack as={Form} id={FORM_ID} noValidate spacing="16px">
          {!isFixed && (
            <CompoundLookup
              workgroupName={workgroupName}
              compounds={compounds}
              onAddNew={setNewCompoundCas}
            />
          )}

          <SelectField
            name="locationId"
            label="Location"
            isRequired
            helpText={
              activeLocations.length === 0
                ? "This workgroup has no locations yet. Add one from Manage locations first."
                : undefined
            }
          >
            <option value="">Choose a location</option>
            {activeLocations.map((l) => (
              <option key={l.id} value={l.id}>
                {l.name}
              </option>
            ))}
          </SelectField>

          <FieldPair>
            <TextField
              name="quantity"
              label="Quantity"
              type="number"
              inputMode="decimal"
              step="any"
              min={0}
              isRequired
            />
            <SelectField name="unit" label="Unit" isRequired>
              <option value="">Choose a unit</option>
              {UNITS.map((u) => (
                <option key={u} value={u}>
                  {u}
                </option>
              ))}
            </SelectField>
          </FieldPair>

          <TextField name="supplier" label="Supplier" />

          <FieldPair>
            <TextField name="lotNumber" label="Lot number" />
            <TextField name="purity" label="Purity" placeholder="e.g. ≥99.5%" />
          </FieldPair>

          <TextField
            name="barcode"
            label="Barcode"
            helpText="Scan or type the label on the container."
          />

          <FieldPair>
            <TextField name="dateReceived" label="Date received" type="date" />
            <TextField name="expiryDate" label="Expiry date" type="date" />
          </FieldPair>
        </VStack>
      </Formik>

      {/*
        Rendered outside the Add stock form on purpose: React events bubble
        through portals, so a submit from inside the form's subtree would also
        submit Add stock.
      */}
      {newCompoundCas && (
        <NewCompoundModal
          workgroupName={workgroupName}
          workbooks={workbooks}
          cas={newCompoundCas}
          onCreated={handleCompoundCreated}
          onClose={() => setNewCompoundCas(null)}
        />
      )}
    </InvModal>
  );
};
