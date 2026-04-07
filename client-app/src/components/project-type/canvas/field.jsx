import { dropTargetForElements } from "@atlaskit/pragmatic-drag-and-drop/element/adapter";
import { Box, Divider, HStack, Text, useToast, VStack } from "@chakra-ui/react";
import { useSectionFields } from "api";
import { Badge } from "components/core/Badge";
import {
  INPUT_TYPES_MAP,
  InputTypePalette,
} from "components/project-type/canvas/field/input-type-palette";
import { INPUT_TYPES_MAP as FIELD_TYPES_MAP } from "components/section-field";
import { STAGES, TOAST_DEFAULTS } from "constants";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useNavigate, useParams, useSearchParams } from "react-router-dom";

import { BASE_PATH } from "./area";
import { FormActions } from "./field/action";
import { FieldManager, Info } from "./field/manager";

export const Field = ({ section, projectType }) => {
  const [searchParams] = useSearchParams();
  const { projectTypeId, sectionTypeId, sectionId } = useParams();
  const navigate = useNavigate();
  const { fields: api } = useBackendApi();

  const canEdit = projectType.stage === STAGES.Draft;
  const isEditing =
    canEdit &&
    searchParams.get("action") === "edit" &&
    searchParams.get("type") === "section-fields";

  const [isLoading, setIsLoading] = useState(false);
  const [isDragOver, setIsDragOver] = useState(false);

  const [fields, setFields] = useState([]);

  const { data, mutate } = useSectionFields(section.id);

  const orderedFields = useMemo(() => {
    return data ? orderFields(data) : [];
  }, [data]);

  useEffect(() => {
    setFields(orderedFields);
  }, [orderedFields]);

  const toast = useToast();
  const dropRef = useRef(null);

  const handleSave = async () => {
    if (fields.length === 0) {
      toast({
        ...TOAST_DEFAULTS,
        title: "No fields to save. Please add fields first.",
        status: "error",
      });
      return;
    }

    try {
      setIsLoading(true);
      const createModel = fields.map(mapFieldToCreateModel);
      await api.save(section.id, createModel);
      toast({
        ...TOAST_DEFAULTS,
        title: "Fields saved",
        status: "success",
      });
      await mutate();
      navigate(
        `${BASE_PATH}/${projectTypeId}/section-types/${sectionTypeId}/sections/${sectionId}`,
        { replace: true },
      );
    } catch (error) {
      console.error(error);
      toast({
        ...TOAST_DEFAULTS,
        status: "error",
        message: error.message,
      });
    } finally {
      setIsLoading(false);
    }
  };

  const handleCancel = () => {
    setFields(orderFields(data));
    navigate(
      `${BASE_PATH}/${projectTypeId}/section-types/${sectionTypeId}/sections/${sectionId}`,
      { replace: true },
    );
  };

  const handleDrop = useCallback(
    ({ source }) => {
      setIsDragOver(false);
      if (!isEditing) {
        toast({
          ...TOAST_DEFAULTS,
          title: "Please enable edit mode first",
          status: "warning",
        });
        return;
      }

      const { inputType } = source.data;

      if (!inputType || !inputType.name) {
        toast({
          ...TOAST_DEFAULTS,
          title: "Invalid input type",
          status: "error",
        });
        return;
      }

      const newField = {
        id: `temp-${Date.now()}`,
        name: `${inputType.name} ${fields.length + 1}`,
        mandatory: false,
        hidden: false,
        inputType,
        sortOrder: fields.length + 1,
        selectFieldOptions: INPUT_TYPES_MAP[inputType.name]?.options || [],
      };

      setFields([...fields, newField]);

      toast({
        ...TOAST_DEFAULTS,
        title: `Added ${inputType.name}`,
        status: "success",
      });
    },
    [fields, toast, isEditing],
  );

  useEffect(() => {
    if (!dropRef.current) return;

    const cleanup = dropTargetForElements({
      element: dropRef.current,
      canDrop: ({ source }) => source.data.type === "input-type",
      onDragEnter: () => isEditing && setIsDragOver(true),
      onDragLeave: () => setIsDragOver(false),
      onDrop: handleDrop,
    });

    return cleanup;
  }, [isEditing, handleDrop]);

  const initialValues = fields.reduce((acc, field) => {
    try {
      acc[field.id] = JSON.parse(field.defaultResponse);
    } catch (e) {
      acc[field.id] = INPUT_TYPES_MAP[field.inputType.name].defaultResponse;
    }
    return acc;
  }, {});

  return (
    <HStack align="start" w="full">
      {isEditing && (
        <InputTypePalette
          onAdd={(inputType) => handleDrop({ source: { data: { inputType } } })}
        />
      )}
      <VStack
        ref={dropRef}
        w="full"
        p={4}
        align="stretch"
        spacing={4}
        borderWidth={1}
        borderRadius={4}
        borderColor={isDragOver ? "green.300" : "gray.200"}
        transition="all 0.2s"
        minH="400px"
        bg={isDragOver ? "green.50" : "white"}
      >
        <HStack justify="space-between" w="full">
          <Badge label="Section fields" colorScheme="orange" />
          {canEdit && (
            <FormActions
              handleSubmit={handleSave}
              handleCancel={handleCancel}
              isEditing={isEditing}
              isLoading={isLoading}
            />
          )}
        </HStack>
        <Divider />

        {fields.length === 0 && <NoFieldsAlert />}

        {Object.keys(initialValues).length > 0 && (
          <Formik enableReinitialize initialValues={initialValues}>
            <Form>
              <VStack spacing={4} align="stretch" w="full">
                {!isEditing ? (
                  fields.map((field) => (
                    <FieldRenderer key={field.id} field={field} />
                  ))
                ) : (
                  <FieldManager fields={fields} setFields={setFields} />
                )}
              </VStack>
            </Form>
          </Formik>
        )}
      </VStack>
    </HStack>
  );
};

const NoFieldsAlert = () => (
  <VStack spacing={4} py={8} align="center" w="full">
    <Box p={6} borderWidth={1} borderStyle="dashed" textAlign="center">
      <Text color="gray.500" fontSize="sm" mb={2}>
        No fields added yet
      </Text>
      <Text color="gray.400" fontSize="xs">
        Enable edit mode to add fields. Drag input types from the palette to add
        them
      </Text>
    </Box>
  </VStack>
);

const FieldRenderer = ({ field, isChild = false, depth = 0 }) => {
  const [, Component] =
    Object.entries(FIELD_TYPES_MAP).find(
      ([key]) => key.toUpperCase() === field.inputType?.name.toUpperCase(),
    ) || [];

  if (!Component) return null;

  const Field = () => (
    <HStack key={field.id} align="start">
      {field.sortOrder != 0 && (
        <Text fontWeight="light" fontSize="xxs" color="gray.500">
          {field.sortOrder}.
        </Text>
      )}
      <VStack spacing={2} align="start" w="full">
        <Component field={field} isDisabled />
        <HStack justify="end" w="full">
          <Info field={field} />
        </HStack>
        <Divider />
      </VStack>
    </HStack>
  );

  return (
    <>
      {isChild ? (
        <Box ml={`${depth * 20}px`} position="relative">
          <Box p={2} borderWidth={1} borderRadius="md" borderColor="purple.200">
            <Field />
          </Box>
        </Box>
      ) : (
        <Field />
      )}

      {field.triggerField && (
        <FieldRenderer field={field.triggerField} depth={depth + 1} isChild />
      )}
    </>
  );
};

/**
 * Maps field to create model
 * @param {*} field - field to map
 * @returns create model
 */
const mapFieldToCreateModel = (field) => {
  const model = {
    id: typeof field.id === "number" ? field.id : null,
    name: field.name,
    inputType: field.inputType.id,
    mandatory: field.mandatory,
    hidden: field.hidden,
    sortOrder: field.sortOrder,
    defaultValue: JSON.stringify(
      INPUT_TYPES_MAP[field.inputType.name].defaultResponse,
    ),
    selectFieldOptions:
      field.selectFieldOptions?.map((option) => ({
        id: typeof option?.id === "number" ? option.id : null,
        name: option.name,
      })) || [],
  };

  if (field.triggerField) {
    model.triggerField = mapFieldToCreateModel(field.triggerField);
    model.triggerValue = field.triggerField.triggerValue;
  }

  return model;
};

/**
 * Orders fields by sortOrder and structures trigger field relationships.
 *
 * @example
 * Input: [
 *   { id: 1, name: "Name", sortOrder: 1 },
 *   { id: 2, name: "Type", sortOrder: 2, triggerField: { id: 3 } },
 *   { id: 3, name: "Details", sortOrder: 0 }
 * ]
 *
 * Output: [
 *   { id: 1, name: "Name", sortOrder: 1 },
 *   { id: 2, name: "Type", sortOrder: 2, triggerField: { id: 3 } }
 * ]
 */
export const orderFields = (fields) => {
  if (!fields?.length) return [];

  const cacheKey = fields
    .map((field) => `${field.id}-${field.triggerField?.id || "none"}`)
    .join(",");

  if (orderFields._cache && orderFields._lastKey === cacheKey)
    return orderFields._cache;

  const fieldMap = new Map(fields.map((field) => [field.id, { ...field }]));
  const childIds = new Set();
  const processed = new Set();

  const populateTriggerField = (field) => {
    if (processed.has(field.id)) return field;

    processed.add(field.id);

    if (field.triggerField?.id && fieldMap.has(field.triggerField.id)) {
      const triggeredField = fieldMap.get(field.triggerField.id);

      field.triggerField = {
        ...populateTriggerField(triggeredField),
        triggerValue: field.triggerField.value,
      };

      childIds.add(triggeredField.id);
    }

    processed.delete(field.id);
    return field;
  };

  const result = [...fieldMap.values()]
    .map(populateTriggerField)
    .filter((field) => !childIds.has(field.id))
    .sort((a, b) => (a.sortOrder || 0) - (b.sortOrder || 0));

  orderFields._cache = result;
  orderFields._lastKey = cacheKey;

  return result;
};

export const DRAG_TYPES = {
  INPUT_TYPE: "input-type",
  FIELD_ITEM: "field-item",
};
