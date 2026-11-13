import { Grid, Text, VStack } from "@chakra-ui/react";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { number, object, string } from "yup";
import { SelectField, TextAreaField, TextField } from "./form-fields.jsx";
import { ErrorMessage, INV, InvModal } from "./ui.jsx";
import { useInventoryAction } from "./use-inventory-action.js";
import { formatDate, formatQuantity } from "./utils.jsx";

/*
 * Modals for the actions on a single container: adjust its quantity, move it
 * to another location, or dispose of it. Each is mounted by the inventory page
 * while open, so its form and error state start fresh every time.
 */

const NOTES_MAX_LENGTH = 500;
const notesSchema = string()
  .trim()
  .max(NOTES_MAX_LENGTH, `Keep notes under ${NOTES_MAX_LENGTH} characters`);

/** Trims optional free text and turns an empty value into null for the API. */
const optionalText = (value) => value.trim() || null;

/**
 * Read-only summary of the container being acted on.
 * @param {object} props
 * @param {object} props.item The InventoryItem.
 * @param {object} props.compound The item's compound (name and CAS).
 */
export const ItemSummary = ({ item, compound }) => (
  <Grid
    as="dl"
    templateColumns="max-content 1fr"
    columnGap="16px"
    rowGap="4px"
    mb="16px"
    px="12px"
    py="10px"
    bg={INV.sub}
    borderRadius="8px"
    fontSize="13px"
    sx={{ dt: { color: INV.muted }, dd: { m: 0 } }}
  >
    <dt>Compound</dt>
    <dd>
      <Text as="b">{compound.name}</Text> ({compound.cas})
    </dd>
    <dt>Location</dt>
    <dd>{item.location.name}</dd>
    <dt>Quantity</dt>
    <dd>
      {formatQuantity(item.quantity)} {item.unit}
    </dd>
    <dt>Lot</dt>
    <dd>{item.lotNumber || "-"}</dd>
    <dt>Expiry</dt>
    <dd>{formatDate(item.expiryDate)}</dd>
  </Grid>
);

/**
 * Set a container's quantity to what was measured, e.g. after a stock check
 * or when some was used outside a recorded reaction.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup the container belongs to.
 * @param {object} props.item The InventoryItem.
 * @param {object} props.compound The item's compound.
 * @param {() => void} props.onClose Closes the modal.
 */
export const AdjustModal = ({ workgroupName, item, compound, onClose }) => {
  const { inventory } = useBackendApi();
  const { run, isLoading, error } = useInventoryAction(workgroupName);

  const schema = object({
    quantity: number()
      .typeError("Enter a number")
      .required("Enter the new quantity")
      .min(0, "Quantity can't be negative")
      .notOneOf([Number(item.quantity)], "Enter a different quantity"),
    notes: notesSchema,
  });

  const handleSubmit = async (values) => {
    const ok = await run(
      () =>
        inventory.adjustItem(workgroupName, item.id, {
          quantity: Number(values.quantity),
          notes: optionalText(values.notes),
        }),
      "Quantity updated",
    );
    if (ok) onClose();
  };

  return (
    <InvModal
      isOpen
      onClose={onClose}
      title="Adjust quantity"
      actionLabel="Save"
      formId="inv-adjust-form"
      isLoading={isLoading}
    >
      {error && <ErrorMessage>{error}</ErrorMessage>}
      <ItemSummary item={item} compound={compound} />
      <Formik
        initialValues={{ quantity: "", notes: "" }}
        validationSchema={schema}
        onSubmit={handleSubmit}
      >
        <VStack as={Form} id="inv-adjust-form" noValidate spacing="16px">
          <TextField
            name="quantity"
            label={`New quantity (${item.unit})`}
            type="number"
            inputMode="decimal"
            step="any"
            min={0}
            helpText="The amount left in the container now."
            isRequired
            autoFocus
          />
          <TextAreaField
            name="notes"
            label="Notes"
            helpText="For example, why the quantity changed."
          />
        </VStack>
      </Formik>
    </InvModal>
  );
};

/**
 * Move a container to another of the workgroup's locations.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup the container belongs to.
 * @param {object} props.item The InventoryItem.
 * @param {object} props.compound The item's compound.
 * @param {object[]} props.locations The workgroup's locations.
 * @param {() => void} props.onClose Closes the modal.
 */
export const MoveModal = ({
  workgroupName,
  item,
  compound,
  locations,
  onClose,
}) => {
  const { inventory } = useBackendApi();
  const { run, isLoading, error } = useInventoryAction(workgroupName);

  const destinations = locations.filter(
    (l) => !l.isArchived && l.id !== item.location.id,
  );

  const schema = object({
    locationId: string().required("Choose a location"),
  });

  const handleSubmit = async (values) => {
    const ok = await run(
      () =>
        inventory.moveItem(workgroupName, item.id, {
          locationId: Number(values.locationId),
        }),
      "Container moved",
    );
    if (ok) onClose();
  };

  // nothing to move to: show an explanation instead of the form
  if (destinations.length === 0) {
    return (
      <InvModal isOpen onClose={onClose} title="Move container">
        <ItemSummary item={item} compound={compound} />
        <Text>
          There are no other locations to move this container to. Add one from
          Manage locations first.
        </Text>
      </InvModal>
    );
  }

  return (
    <InvModal
      isOpen
      onClose={onClose}
      title="Move container"
      actionLabel="Move"
      formId="inv-move-form"
      isLoading={isLoading}
    >
      {error && <ErrorMessage>{error}</ErrorMessage>}
      <ItemSummary item={item} compound={compound} />
      <Formik
        initialValues={{ locationId: "" }}
        validationSchema={schema}
        onSubmit={handleSubmit}
      >
        <Form id="inv-move-form" noValidate>
          <SelectField name="locationId" label="New location" isRequired>
            <option value="">Choose a location</option>
            {destinations.map((l) => (
              <option key={l.id} value={l.id}>
                {l.name}
              </option>
            ))}
          </SelectField>
        </Form>
      </Formik>
    </InvModal>
  );
};

/**
 * Confirm disposing of a container. Disposal empties and archives it; its
 * history is kept.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup the container belongs to.
 * @param {object} props.item The InventoryItem.
 * @param {object} props.compound The item's compound.
 * @param {() => void} props.onClose Closes the modal.
 */
export const DisposeModal = ({ workgroupName, item, compound, onClose }) => {
  const { inventory } = useBackendApi();
  const { run, isLoading, error } = useInventoryAction(workgroupName);

  const handleSubmit = async (values) => {
    const ok = await run(
      () =>
        inventory.disposeItem(workgroupName, item.id, {
          notes: optionalText(values.notes),
        }),
      "Container disposed",
    );
    if (ok) onClose();
  };

  return (
    <InvModal
      isOpen
      onClose={onClose}
      title="Dispose of container"
      actionLabel="Dispose"
      formId="inv-dispose-form"
      isLoading={isLoading}
      danger
    >
      {error && <ErrorMessage>{error}</ErrorMessage>}
      <ItemSummary item={item} compound={compound} />
      <Text mb="16px">
        The container's quantity will be set to zero and it will be archived.
        Its history is kept, and you can still see it with Show archived.
      </Text>
      <Formik
        initialValues={{ notes: "" }}
        validationSchema={object({ notes: notesSchema })}
        onSubmit={handleSubmit}
      >
        <Form id="inv-dispose-form" noValidate>
          <TextAreaField
            name="notes"
            label="Reason"
            helpText="For example, expired, contaminated or empty."
          />
        </Form>
      </Formik>
    </InvModal>
  );
};
