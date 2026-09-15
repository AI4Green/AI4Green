import {
  draggable,
  dropTargetForElements,
} from "@atlaskit/pragmatic-drag-and-drop/element/adapter";
import { reorder } from "@atlaskit/pragmatic-drag-and-drop/reorder";
import {
  Button,
  FormControl,
  FormErrorMessage,
  HStack,
  IconButton,
  Input,
  Stack,
  Text,
  VStack,
} from "@chakra-ui/react";
import { useField } from "formik";
import { useEffect, useRef, useState } from "react";
import { FaPlus, FaRegTimesCircle } from "react-icons/fa";
import { RxDragHandleDots2 } from "react-icons/rx";

import { generateUniqueId, useDropMonitor } from "./draggable-list-field";

const DRAG_TYPE = "draggable-list-item";

export const InlineDraggableListField = ({
  label,
  name,
  isDisabled,
  addLabel = "Add new item",
}) => {
  const [field, , helpers] = useField(name);

  const handleAdd = () => {
    const newId = generateUniqueId(field.value.map((i) => i.id));
    helpers.setValue([
      ...field.value,
      { id: newId, order: field.value.length + 1, content: "" },
    ]);
  };

  const handleDelete = (id) => {
    helpers.setValue(
      field.value
        .filter((i) => i.id !== id)
        .map((i, idx) => ({ ...i, order: idx + 1 })),
    );
  };

  useDropMonitor(isDisabled, (from, to) => {
    if (from === to) return;
    const reordered = reorder({
      list: [...field.value].sort((a, b) => a.order - b.order),
      startIndex: from,
      finishIndex: to,
    }).map((item, idx) => ({ ...item, order: idx + 1 }));
    helpers.setValue(reordered);
  });

  return (
    <Stack align="start" spacing={4} w="100%">
      <Text as="b">{label}</Text>
      <VStack w="full" align="start" spacing={4}>
        {field.value
          .sort((a, b) => a.order - b.order)
          .map((item, idx) => (
            <Item
              key={item.id}
              item={item}
              index={idx}
              name={name}
              onDelete={handleDelete}
              isDisabled={isDisabled}
            />
          ))}
      </VStack>

      <Button
        size="xs"
        variant="outline"
        colorScheme="green"
        leftIcon={<FaPlus />}
        onClick={handleAdd}
        isDisabled={isDisabled}
      >
        {addLabel}
      </Button>
    </Stack>
  );
};

const Item = ({ item, index, name, onDelete, isDisabled }) => {
  const [isDragging, setIsDragging] = useState(false);
  const [isOver, setIsOver] = useState(false);
  const dragRef = useRef(null);
  const dropRef = useRef(null);

  useEffect(() => {
    if (isDisabled) return;

    const cleanupDrag = draggable({
      element: dragRef.current,
      getInitialData: () => ({ itemId: item.id, index, type: DRAG_TYPE }),
      onDragStart: () => setIsDragging(true),
      onDrop: () => setIsDragging(false),
    });

    const cleanupDrop = dropTargetForElements({
      element: dropRef.current,
      getData: () => ({ index, type: DRAG_TYPE }),
      canDrop: ({ source }) => source.data.type === DRAG_TYPE,
      onDragEnter: () => setIsOver(true),
      onDragLeave: () => setIsOver(false),
      onDrop: () => setIsOver(false),
    });

    return () => {
      cleanupDrag();
      cleanupDrop();
    };
  }, [index, isDisabled, item.id]);

  return (
    <HStack
      ref={dropRef}
      w="100%"
      p={2}
      border="1px solid"
      borderColor={isOver ? "blue.300" : "gray.200"}
      borderRadius="md"
      align="center"
      justify="space-between"
      gap={2}
      bg={isDragging ? "blue.50" : "white"}
      boxShadow={isDragging ? "md" : "sm"}
      opacity={isDragging ? 0.7 : 1}
      transition="all 0.2s"
    >
      <HStack ref={dragRef} cursor="grab" _active={{ cursor: "grabbing" }}>
        {!isDisabled && <RxDragHandleDots2 />}
        <Text
          fontWeight="light"
          mr={2}
          fontSize="xs"
          color={isDisabled ? "gray.500" : "gray.800"}
        >
          {item.order}.
        </Text>
      </HStack>

      <InputField name={`${name}[${index}].content`} isDisabled={isDisabled} />

      <IconButton
        aria-label="Remove item"
        icon={<FaRegTimesCircle />}
        size="sm"
        variant="ghost"
        colorScheme="red"
        onClick={() => onDelete(item.id)}
        isDisabled={isDisabled}
      />
    </HStack>
  );
};

const InputField = ({ name, ...p }) => {
  const [field, meta] = useField(name);

  return (
    <FormControl id={field.name} isInvalid={meta.error && meta.touched}>
      <Input {...field} {...p} />
      {meta.error && meta.touched && (
        <FormErrorMessage>
          <Text fontSize="xs">{meta.error}</Text>
        </FormErrorMessage>
      )}
    </FormControl>
  );
};
