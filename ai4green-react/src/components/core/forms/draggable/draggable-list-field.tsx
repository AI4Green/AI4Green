import {
  draggable,
  dropTargetForElements,
  monitorForElements,
} from "@atlaskit/pragmatic-drag-and-drop/element/adapter";
import { reorder } from "@atlaskit/pragmatic-drag-and-drop/reorder";
import {
  Box,
  Button,
  Flex,
  IconButton,
  Stack,
  Text,
  useDisclosure,
  VStack,
} from "@chakra-ui/react";
import { Modal } from "components/core/modal.tsx";
import { Field, Form, Formik, useField } from "formik";
import { useEffect, useRef, useState } from "react";
import { FaBook, FaEdit, FaRegTimesCircle } from "react-icons/fa";
import { object, string } from "yup";

import { TextAreaField } from "..";

const DRAG_TYPE = "draggable-list-item";

const validationSchema = () =>
  object().shape({
    itemInput: string().required("Input required"),
  });

const ItemModal = ({ isOpen, onClose, onSubmit, title, initialContent }) => {
  const formRef = useRef();
  return (
    <Modal
      title={title}
      isOpen={isOpen}
      onClose={onClose}
      onAction={() => formRef.current?.handleSubmit()}
      body={
        <Formik
          enableReinitialize
          innerRef={formRef}
          initialValues={{ itemInput: initialContent || "" }}
          validationSchema={validationSchema}
          onSubmit={(values) => {
            onSubmit(values.itemInput);
            onClose();
          }}
        >
          <Form noValidate>
            <TextAreaField name="itemInput" placeholder="Enter here" />
          </Form>
        </Formik>
      }
    />
  );
};

const Item = ({ item, index, onEdit, onDelete, isDisabled }) => {
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
    <Box
      ref={dropRef}
      w="100%"
      borderRadius="md"
      border="2px solid"
      borderColor={isOver ? "blue.300" : "transparent"}
      bg={isOver ? "blue.50" : "transparent"}
    >
      <Flex
        ref={dragRef}
        px={2}
        py={2}
        my={1}
        borderBottomWidth={1}
        align="center"
        opacity={isDragging ? 0.5 : 1}
        cursor={isDisabled ? "default" : "grab"}
        _active={{ cursor: isDisabled ? "default" : "grabbing" }}
      >
        <Text fontWeight="semibold" mr={2} fontSize="xs">
          {item.order}.
        </Text>
        <Text fontSize="sm" flex="1">
          {item.content}
        </Text>
        <IconButton
          aria-label="Edit item"
          icon={<FaEdit />}
          size="sm"
          variant="ghost"
          onClick={() => onEdit(item)}
          isDisabled={isDisabled}
        />
        <IconButton
          aria-label="Remove item"
          icon={<FaRegTimesCircle />}
          size="sm"
          variant="ghost"
          colorScheme="red"
          onClick={() => onDelete(item.id)}
          isDisabled={isDisabled}
        />
      </Flex>
    </Box>
  );
};

export const DraggableListField = ({
  label,
  name,
  isDisabled,
  addLabel = "Add new item",
  editLabel = "Edit Item",
}) => {
  const [field, , helpers] = useField(name);
  const [editItem, setEditItem] = useState(null);
  const addModal = useDisclosure();
  const editModal = useDisclosure();

  const handleAdd = (content) => {
    const newId = generateUniqueId(field.value.map((i) => i.id));
    helpers.setValue([
      ...field.value,
      { id: newId, order: field.value.length + 1, content },
    ]);
  };

  const handleEdit = (content) => {
    helpers.setValue(
      field.value.map((i) => (i.id === editItem.id ? { ...i, content } : i)),
    );
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
    <Field name={name}>
      {() => (
        <Stack align="start" spacing={4} w="100%" py={2}>
          <Text as="b">{label}</Text>

          {field.value?.length > 0 && (
            <VStack
              w="100%"
              align="start"
              spacing={2}
              p={4}
              borderWidth={1}
              borderRadius="md"
            >
              {field.value
                .sort((a, b) => a.order - b.order)
                .map((item, idx) => (
                  <Item
                    key={item.id}
                    item={item}
                    index={idx}
                    isDisabled={isDisabled}
                    onEdit={(i) => {
                      setEditItem(i);
                      editModal.onOpen();
                    }}
                    onDelete={handleDelete}
                  />
                ))}
            </VStack>
          )}

          <Button
            size="sm"
            variant="outline"
            colorScheme="blue"
            leftIcon={<FaBook />}
            onClick={addModal.onOpen}
            isDisabled={isDisabled}
          >
            {addLabel}
          </Button>

          <ItemModal
            title={`📖 ${addLabel}`}
            isOpen={addModal.isOpen}
            onClose={addModal.onClose}
            onSubmit={handleAdd}
          />

          {editItem && (
            <ItemModal
              title={`📖 ${editLabel}`}
              isOpen={editModal.isOpen}
              onClose={editModal.onClose}
              onSubmit={handleEdit}
              initialContent={editItem.content}
            />
          )}
        </Stack>
      )}
    </Field>
  );
};

export const generateUniqueId = (ids) => {
  let id;
  do {
    id = `temp-${Math.random().toString(36).slice(2, 11)}`;
  } while (ids.includes(id));
  return id;
};

export const useDropMonitor = (isDisabled, onReorder) => {
  useEffect(() => {
    if (isDisabled) return;
    return monitorForElements({
      onDrop: ({ source, location }) => {
        const dest = location.current.dropTargets[0];
        if (!dest) return;
        const { index: from } = source.data;
        const { index: to } = dest.data;
        if (typeof from === "number" && typeof to === "number") {
          onReorder(from, to);
        }
      },
    });
  }, [isDisabled, onReorder]);
};
