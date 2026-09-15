import {
  Box,
  HStack,
  IconButton,
  Text,
  useDisclosure,
  VStack,
} from "@chakra-ui/react";
import { useInputTypes } from "api/field";
import { ActionButton } from "components/core/action-button";
import { FormikInput, MultiSelectField, Switch } from "components/core/forms";
import { Modal } from "components/core/modal";
import { INPUT_TYPES_MAP as FIELD_TYPES_MAP } from "components/section-field";
import { INPUT_TYPES } from "constants";
import { Form, Formik } from "formik";
import { capitalise } from "helpers/strings";
import { useRef, useState } from "react";
import { FaPencilAlt, FaPlus, FaSave, FaTrash } from "react-icons/fa";
import { TbCancel } from "react-icons/tb";
import { useSearchParams } from "react-router-dom";
import CreatableSelect from "react-select/creatable";
import { array, object, string } from "yup";

import { INPUT_TYPES_MAP } from "./input-type-palette";

export const FieldActions = ({ field, fields, setFields, isChild }) => {
  const {
    isOpen: isOpenEdit,
    onOpen: onOpenEdit,
    onClose: onCloseEdit,
  } = useDisclosure();

  const {
    isOpen: isOpenAddChild,
    onOpen: onOpenAddChild,
    onClose: onCloseAddChild,
  } = useDisclosure();

  const handleDelete = () => {
    const remove = (f) => {
      if (f.id === field.id) return null;

      if (f.triggerField) return { ...f, triggerField: remove(f.triggerField) };

      return f;
    };

    const model = fields.map(remove).filter((f) => f !== null);
    setFields(model);
  };

  const handleEditSubmit = (values) => {
    const update = (f) => {
      if (f.id === field.id) {
        return {
          ...f,
          ...values,
          sortOrder: values.hidden ? 0 : f.sortOrder,
        };
      }

      if (f.triggerField) return { ...f, triggerField: update(f.triggerField) };

      return f;
    };

    const model = fields.map(update);
    setFields(model);
    onCloseEdit();
  };

  const handleAddChildSubmit = (values) => {
    const { inputType, triggerValue } = values;

    const child = {
      id: `temp-child-${Date.now()}`,
      name: `${field.name} child - ${inputType.name}`,
      mandatory: false,
      hidden: true,
      inputType,
      sortOrder: 0,
      selectFieldOptions: INPUT_TYPES_MAP[inputType.name]?.options || [],
      triggerValue,
    };

    const add = (f) => {
      if (f.id === field.id) return { ...f, triggerField: child };

      if (f.triggerField) return { ...f, triggerField: add(f.triggerField) };

      return f;
    };

    const model = fields.map(add);
    setFields(model);
    onCloseAddChild();
  };

  const actions = {
    delete: {
      isEligible: () => true,
      icon: <FaTrash />,
      label: "Delete",
      onClick: handleDelete,
    },
    edit: {
      isEligible: () => true,
      icon: <FaPencilAlt />,
      label: "Edit",
      onClick: onOpenEdit,
    },
    addChild: {
      isEligible: () => !field.triggerField,
      icon: <FaPlus />,
      label: "Add Child",
      onClick: onOpenAddChild,
    },
  };

  return (
    <HStack justify="end">
      <ActionButton size="xs" actions={actions} />
      {isOpenEdit && (
        <FieldEditModal
          isOpen={isOpenEdit}
          onClose={onCloseEdit}
          field={field}
          handleEditSubmit={handleEditSubmit}
          isChild={isChild}
        />
      )}
      {isOpenAddChild && (
        <FieldAddChildModal
          isOpen={isOpenAddChild}
          onClose={onCloseAddChild}
          field={field}
          handleAddChildSubmit={handleAddChildSubmit}
        />
      )}
    </HStack>
  );
};

const FieldEditModal = ({
  isOpen,
  onClose,
  field,
  handleEditSubmit,
  isChild,
}) => {
  const selectFieldOptions = field.selectFieldOptions?.map((option) => ({
    id: option.id,
    label: option.name,
    value: option.name,
  }));
  const [selectedOptions, setSelectedOptions] = useState(selectFieldOptions);

  const initialValues = {
    name: field.name,
    mandatory: field.mandatory,
    hidden: field.hidden,
    triggerValue: isChild ? field.triggerValue : null,
  };

  const validationSchema = object({
    name: string().required("Field name is required."),
    ...(isChild && {
      triggerValue: string().required("Trigger value is required."),
    }),
  });

  const formRef = useRef();
  const modalBody = (
    <Formik
      enableReinitialize
      innerRef={formRef}
      initialValues={initialValues}
      onSubmit={(values, { resetForm }) => {
        values.selectFieldOptions = selectedOptions?.map((option) => ({
          id: option.id || `new-${option.value}`,
          name: option.value,
        }));
        handleEditSubmit(values);
        resetForm();
      }}
      validationSchema={validationSchema}
    >
      <Form noValidate>
        <VStack spacing={8} align="stretch">
          <FormikInput name="name" label="Name" isRequired />
          <HStack>
            <Switch name="mandatory" label="Mandatory" colorScheme="orange" />
            {!isChild && (
              <Switch name="hidden" label="Hidden" colorScheme="blue" />
            )}
          </HStack>

          {field.inputType.name === INPUT_TYPES.Multiple ||
          field.inputType.name === INPUT_TYPES.Radio ? (
            <HStack>
              <Text fontSize="sm">Options</Text>
              <Box flex={1}>
                <CreatableSelect
                  isCreatable
                  isMulti
                  options={selectFieldOptions}
                  value={selectedOptions}
                  onChange={(value) => {
                    setSelectedOptions(value);
                  }}
                  placeholder="Type to add new option"
                />
              </Box>
            </HStack>
          ) : null}

          {isChild && (
            <FormikInput name="triggerValue" label="Trigger Value" isRequired />
          )}
        </VStack>
      </Form>
    </Formik>
  );

  return (
    <Modal
      isOpen={isOpen}
      onClose={onClose}
      title="Edit"
      body={modalBody}
      onAction={() => formRef.current.handleSubmit()}
      actionBtnCaption="Update"
    />
  );
};

const FieldAddChildModal = ({ isOpen, onClose, handleAddChildSubmit }) => {
  const { data: inputTypes } = useInputTypes();

  const validationSchema = object({
    inputType: array()
      .of(string())
      .min(1, "Input type is required.")
      .required("Input type is required."),
    triggerValue: string().required("Trigger value is required."),
  });

  const formRef = useRef();
  const modalBody = (
    <Formik
      enableReinitialize
      innerRef={formRef}
      initialValues={{
        inputType: [],
        triggerValue: "",
      }}
      validationSchema={validationSchema}
      onSubmit={(values, { resetForm }) => {
        const selectedInputType =
          values.inputType.length > 0
            ? inputTypes.find(
                (inputType) => String(inputType.id) === values.inputType[0],
              )
            : null;

        if (selectedInputType) {
          handleAddChildSubmit({
            inputType: selectedInputType,
            triggerValue: values.triggerValue,
          });
        }
        resetForm();
      }}
    >
      {({ values }) => {
        const selectedInputType = inputTypes?.find(
          (inputType) => String(inputType.id) === values.inputType[0],
        );
        const [, Component] =
          Object.entries(FIELD_TYPES_MAP).find(
            ([key]) =>
              key.toUpperCase() === selectedInputType?.name.toUpperCase(),
          ) || [];

        return (
          <Form noValidate>
            <VStack spacing={8} align="stretch">
              <FormikInput
                name="triggerValue"
                label="Trigger Value"
                isRequired
              />

              <VStack align="stretch">
                <MultiSelectField
                  name="inputType"
                  label="Input Type"
                  isRequired
                  options={inputTypes?.map((inputType) => ({
                    label:
                      INPUT_TYPES_MAP[inputType.name]?.label ||
                      capitalise(inputType.name),
                    value: String(inputType.id),
                  }))}
                />

                {Component && (
                  <>
                    <Component
                      field={{
                        id: "field-preview",
                        name: "",
                        selectFieldOptions:
                          INPUT_TYPES_MAP[selectedInputType.name]?.options ||
                          [],
                        inputType: selectedInputType,
                      }}
                      isDisabled
                    />
                    <Text fontSize="xxs" fontWeight="thin" color="gray.500">
                      Preview may not always be available due to the nature of
                      the input type.
                    </Text>
                  </>
                )}
              </VStack>
            </VStack>
          </Form>
        );
      }}
    </Formik>
  );
  return (
    <Modal
      isOpen={isOpen}
      onClose={onClose}
      title="Add Child"
      body={modalBody}
      onAction={() => formRef.current.handleSubmit()}
      actionBtnCaption="Add"
    />
  );
};

export const FormActions = ({
  isEditing,
  isLoading,
  handleSubmit,
  handleCancel,
}) => {
  const [, setSearchParams] = useSearchParams();

  if (isEditing) {
    return (
      <HStack spacing={4}>
        <IconButton
          size="sm"
          icon={<FaSave />}
          aria-label="Save"
          variant="ghost"
          colorScheme="green"
          onClick={handleSubmit}
          fontSize="lg"
          isLoading={isLoading}
        />
        <IconButton
          size="sm"
          fontSize="lg"
          icon={<TbCancel />}
          aria-label="Cancel"
          variant="ghost"
          colorScheme="yellow"
          onClick={handleCancel}
          isLoading={isLoading}
        />
      </HStack>
    );
  }

  return (
    <IconButton
      size="xs"
      icon={<FaPencilAlt />}
      aria-label="Edit"
      variant="ghost"
      colorScheme="blue"
      onClick={() => {
        setSearchParams({ action: "edit", type: "section-fields" });
      }}
    />
  );
};
