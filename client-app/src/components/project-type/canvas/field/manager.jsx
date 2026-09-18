import {
  draggable,
  dropTargetForElements,
  monitorForElements,
} from "@atlaskit/pragmatic-drag-and-drop/element/adapter";
import { reorder } from "@atlaskit/pragmatic-drag-and-drop/reorder";
import { Box, HStack, Icon, Text, VStack } from "@chakra-ui/react";
import { Badge } from "components/core/Badge";
import { DRAG_TYPES } from "components/project-type/canvas/field";
import { INPUT_TYPES_MAP as FIELD_TYPES_MAP } from "components/section-field";
import { useCallback, useEffect, useRef, useState } from "react";
import {
  FaArrowRight,
  FaCheckCircle,
  FaEyeSlash,
  FaGripVertical,
} from "react-icons/fa";

import { FieldActions } from "./action";
import { INPUT_TYPES_MAP } from "./input-type-palette";

export const FieldManager = ({ fields, setFields }) => {
  const handleReorder = useCallback(
    (from, to) => {
      if (from === to) return;
      const reordered = reorder({
        list: [...fields].sort((a, b) => a.sortOrder - b.sortOrder),
        startIndex: from,
        finishIndex: to,
      }).map((item, idx) => ({ ...item, sortOrder: idx + 1 }));
      setFields(reordered);
    },
    [fields, setFields],
  );

  useEffect(() => {
    return monitorForElements({
      onDrop: ({ source, location }) => {
        const dest = location.current.dropTargets[0];
        if (!dest) return;
        handleReorder(source.data.index, dest.data.index);
      },
    });
  }, [handleReorder]);

  return (
    <>
      {fields.map((field, index) => (
        <FieldItem
          key={field.id}
          field={field}
          index={index}
          fields={fields}
          setFields={setFields}
        />
      ))}
    </>
  );
};

const FieldItem = ({ field, index, fields, setFields, depth = 0 }) => {
  const isChild = isChildField(fields, field.id);
  const dragAndDropProps = useDragAndDrop(field, index, isChild);

  return (
    <>
      {isChild ? (
        <ChildField
          field={field}
          fields={fields}
          setFields={setFields}
          depth={depth}
        />
      ) : (
        <ParentField
          field={field}
          fields={fields}
          setFields={setFields}
          dragAndDropProps={dragAndDropProps}
        />
      )}

      {field.triggerField && (
        <FieldItem
          field={field.triggerField}
          index={-1}
          fields={fields}
          setFields={setFields}
          depth={depth + 1}
        />
      )}
    </>
  );
};

const ParentField = ({ field, fields, setFields, dragAndDropProps }) => {
  const { isDragging, isOver, dropRef } = dragAndDropProps;

  return (
    <VStack
      ref={dropRef}
      p={3}
      borderWidth={1}
      borderRadius="md"
      borderColor={isOver ? "blue.300" : "gray.200"}
      bg={isOver ? "blue.50" : isDragging ? "gray.100" : "white"}
      align="stretch"
      transition="all 0.1s"
      w="full"
    >
      <FieldContent
        field={field}
        fields={fields}
        setFields={setFields}
        dragRef={dragAndDropProps.dragRef}
        isChild={false}
      />
    </VStack>
  );
};

const ChildField = ({ field, fields, setFields, depth }) => (
  <Box ml={`${depth * 20}px`} position="relative">
    <VStack
      p={3}
      borderWidth={1}
      borderRadius="md"
      borderColor="purple.200"
      bg="purple.25"
      align="stretch"
      w="full"
    >
      <FieldContent
        field={field}
        fields={fields}
        setFields={setFields}
        isChild
      />
    </VStack>
  </Box>
);

const FieldContent = ({ field, fields, setFields, dragRef, isChild }) => {
  const [, Component] =
    Object.entries(FIELD_TYPES_MAP).find(
      ([key]) => key.toUpperCase() === field.inputType?.title.toUpperCase(),
    ) || [];

  if (!Component) return null;

  return (
    <>
      <HStack justify="space-between">
        <Text fontWeight="light" fontSize="xs" color="gray.500">
          {field.sortOrder != 0 && `${field.sortOrder}.`}
        </Text>
        <Info field={field} />
      </HStack>
      <HStack>
        {!isChild && (
          <VStack ref={dragRef} cursor="grab" _active={{ cursor: "grabbing" }}>
            <Icon as={FaGripVertical} color="gray.400" fontSize="xl" />
          </VStack>
        )}

        <Component field={field} isDisabled />
      </HStack>
      <FieldActions
        field={field}
        fields={fields}
        setFields={setFields}
        isChild={isChild}
      />
    </>
  );
};

export const Info = ({ field }) => (
  <HStack>
    {field.triggerValue && (
      <Badge
        label={`Trigger Cause: ${field.triggerValue}`}
        leftIcon={FaArrowRight}
        colorScheme="purple"
        variant="outline"
        fontSize="xxs"
        fontWeight="light"
      />
    )}
    {field.mandatory && (
      <Badge
        label="Mandatory"
        leftIcon={FaCheckCircle}
        colorScheme="orange"
        variant="outline"
        fontSize="xxs"
        fontWeight="light"
      />
    )}
    {field.hidden && (
      <Badge
        label="Hidden"
        leftIcon={FaEyeSlash}
        colorScheme="blue"
        variant="outline"
        fontSize="xxs"
        fontWeight="light"
      />
    )}

    <Badge
      label={INPUT_TYPES_MAP[field.inputType.title].label}
      leftIcon={INPUT_TYPES_MAP[field.inputType.title].icon}
      colorScheme="gray"
      fontSize="xxs"
      fontWeight="light"
    />
  </HStack>
);

const useDragAndDrop = (field, index, isChild) => {
  const [isDragging, setIsDragging] = useState(false);
  const [isOver, setIsOver] = useState(false);
  const dragRef = useRef(null);
  const dropRef = useRef(null);

  useEffect(() => {
    if (isChild || !dragRef.current || !dropRef.current) return;

    const cleanupDrag = draggable({
      element: dragRef.current,
      getInitialData: () => ({
        itemId: field.id,
        index,
        type: DRAG_TYPES.FIELD_ITEM,
      }),
      onDragStart: () => setIsDragging(true),
      onDrop: () => setIsDragging(false),
    });

    const cleanupDrop = dropTargetForElements({
      element: dropRef.current,
      getData: () => ({ index, type: DRAG_TYPES.FIELD_ITEM }),
      canDrop: ({ source }) => source.data.type === DRAG_TYPES.FIELD_ITEM,
      onDragEnter: () => setIsOver(true),
      onDragLeave: () => setIsOver(false),
      onDrop: () => setIsOver(false),
    });

    return () => {
      cleanupDrag();
      cleanupDrop();
    };
  }, [index, field.id, isChild]);

  return { isDragging, isOver, dragRef, dropRef };
};

/**
 * Checks if field is child at any depth in the collection
 * @param {*} fields - fields collection
 * @param {*} targetId - target field id
 * @returns boolean
 */
const isChildField = (fields, targetId) => {
  for (const f of fields) {
    if (f.triggerField?.id === targetId) {
      return true;
    }
    if (f.triggerField && isChildField([f.triggerField], targetId)) {
      return true;
    }
  }
  return false;
};
