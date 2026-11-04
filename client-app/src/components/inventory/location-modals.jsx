import { Text, VStack } from "@chakra-ui/react";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { object, string } from "yup";
import { TextAreaField, TextField } from "./form-fields.jsx";
import { ErrorMessage, InvModal } from "./ui.jsx";
import { useInventoryAction } from "./use-inventory-action.js";

const FORM_ID = "inv-location-form";

/**
 * Builds the validation schema. Location names must be unique within the
 * workgroup, ignoring case.
 * @param {object[]} locations The workgroup's locations.
 * @param {object} [current] The location being edited, excluded from the check.
 */
const buildSchema = (locations, current) => {
  const takenNames = locations
    .filter((l) => l.id !== current?.id)
    .map((l) => l.name.trim().toLowerCase());

  return object({
    name: string()
      .trim()
      .required("Enter a name")
      .max(100, "Keep the name under 100 characters")
      .test(
        "unique-name",
        "There's already a location with this name",
        (value) => !takenNames.includes((value ?? "").trim().toLowerCase()),
      ),
    description: string()
      .trim()
      .max(255, "Keep the description under 255 characters"),
  });
};

/**
 * Add a storage location, or edit one when `location` is passed.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup the location belongs to.
 * @param {object[]} props.locations The workgroup's locations.
 * @param {object} [props.location] The location to edit.
 * @param {() => void} props.onClose Closes the modal.
 */
export const LocationFormModal = ({
  workgroupName,
  locations,
  location,
  onClose,
}) => {
  const { inventory } = useBackendApi();
  const { run, isLoading, error } = useInventoryAction(workgroupName);
  const isEdit = Boolean(location);

  const handleSubmit = async (values) => {
    const payload = {
      name: values.name.trim(),
      description: values.description.trim() || null,
    };
    const ok = await run(
      () =>
        isEdit
          ? inventory.editLocation(workgroupName, location.id, payload)
          : inventory.addLocation(workgroupName, payload),
      isEdit ? "Location updated" : "Location added",
    );
    if (ok) onClose();
  };

  return (
    <InvModal
      isOpen
      onClose={onClose}
      title={isEdit ? "Edit location" : "Add location"}
      actionLabel={isEdit ? "Save" : "Add"}
      formId={FORM_ID}
      isLoading={isLoading}
    >
      {error && <ErrorMessage>{error}</ErrorMessage>}
      <Formik
        initialValues={{
          name: location?.name ?? "",
          description: location?.description ?? "",
        }}
        validationSchema={buildSchema(locations, location)}
        onSubmit={handleSubmit}
      >
        <VStack as={Form} id={FORM_ID} noValidate spacing="16px">
          <TextField
            name="name"
            label="Name"
            placeholder="e.g. Flammables cabinet A"
            isRequired
            autoFocus
          />
          <TextAreaField
            name="description"
            label="Description"
            helpText="Where to find it, e.g. the room and which cabinet or shelf."
          />
        </VStack>
      </Formik>
    </InvModal>
  );
};

/**
 * Confirm archiving or restoring a location. A location that still holds
 * containers can't be archived, so the modal explains that instead.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup the location belongs to.
 * @param {object} props.location The location to archive or restore.
 * @param {number} props.containerCount Active containers in the location.
 * @param {() => void} props.onClose Closes the modal.
 */
export const ArchiveLocationModal = ({
  workgroupName,
  location,
  containerCount,
  onClose,
}) => {
  const { inventory } = useBackendApi();
  const { run, isLoading, error } = useInventoryAction(workgroupName);
  const isRestore = location.isArchived;

  const handleConfirm = async () => {
    const ok = await run(
      () =>
        inventory.setLocationArchived(workgroupName, location.id, !isRestore),
      isRestore ? "Location restored" : "Location archived",
    );
    if (ok) onClose();
  };

  if (!isRestore && containerCount > 0) {
    return (
      <InvModal isOpen onClose={onClose} title="Archive location">
        <Text>
          <Text as="b">{location.name}</Text> still holds {containerCount}{" "}
          {containerCount === 1 ? "container" : "containers"}. Move or dispose
          of {containerCount === 1 ? "it" : "them"} before archiving the
          location.
        </Text>
      </InvModal>
    );
  }

  return (
    <InvModal
      isOpen
      onClose={onClose}
      title={isRestore ? "Restore location" : "Archive location"}
      actionLabel={isRestore ? "Restore" : "Archive"}
      onAction={handleConfirm}
      isLoading={isLoading}
      danger={!isRestore}
    >
      {error && <ErrorMessage>{error}</ErrorMessage>}
      <Text>
        {isRestore ? (
          <>
            Restore <Text as="b">{location.name}</Text>? Containers can be
            stored there again.
          </>
        ) : (
          <>
            Archive <Text as="b">{location.name}</Text>? It will no longer be
            offered when adding or moving containers. You can restore it later.
          </>
        )}
      </Text>
    </InvModal>
  );
};
