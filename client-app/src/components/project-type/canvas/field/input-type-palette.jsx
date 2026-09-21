import { draggable } from "@atlaskit/pragmatic-drag-and-drop/element/adapter";
import { Divider, HStack, Icon, Text, VStack } from "@chakra-ui/react";
import { useInputTypes } from "api";
import { Badge } from "components/core/Badge";
import { DRAG_TYPES } from "components/project-type/canvas/field";
import { INPUT_TYPES } from "constants";
import { useEffect, useRef, useState } from "react";
import {
  FaAlignLeft,
  FaCalendarAlt,
  FaChartBar,
  FaChartLine,
  FaCheckCircle,
  FaEdit,
  FaExclamationTriangle,
  FaFile,
  FaFileUpload,
  FaFlask,
  FaHashtag,
  FaHeading,
  FaListUl,
  FaSort,
  FaTable,
  FaTrash,
} from "react-icons/fa";

export const INPUT_TYPES_MAP = {
  [INPUT_TYPES.Content]: {
    label: "Content",
    icon: FaAlignLeft,
    description: "Static text content",
    defaultResponse: "",
  },
  [INPUT_TYPES.Text]: {
    label: "Text",
    icon: FaFile,
    description: "Text input field",
    defaultResponse: "",
  },
  [INPUT_TYPES.Description]: {
    label: "Description",
    icon: FaAlignLeft,
    description: "Multi-line text area",
    defaultResponse: "",
  },
  [INPUT_TYPES.Number]: {
    label: "Number",
    icon: FaHashtag,
    description: "Numeric input field",
    defaultResponse: 0,
  },
  [INPUT_TYPES.ReactionScheme]: {
    label: "Reaction Scheme",
    icon: FaFlask,
    description: "Reaction scheme input field",
    defaultResponse: [],
  },
  [INPUT_TYPES.MultiReactionScheme]: {
    label: "Multi Reaction Scheme",
    icon: FaFlask,
    description: "Multi reaction scheme input field",
    defaultResponse: [],
  },
  [INPUT_TYPES.Radio]: {
    label: "Radio",
    icon: FaCheckCircle,
    description: "Radio button field",
    options: [
      { name: "No", value: "No" },
      { name: "Yes", value: "Yes" },
    ],
    defaultResponse: [],
  },
  [INPUT_TYPES.Multiple]: {
    label: "Multiple",
    icon: FaListUl,
    description: "Multiple choice field",
    defaultResponse: [],
  },
  [INPUT_TYPES.File]: {
    label: "File",
    icon: FaFileUpload,
    description: "File upload field",
    defaultResponse: [],
  },
  [INPUT_TYPES.Header]: {
    label: "Header",
    icon: FaHeading,
    description: "Static header text",
    defaultResponse: "",
  },
  [INPUT_TYPES.ChemicalDisposalTable]: {
    label: "Chemical Disposal Table",
    icon: FaTrash,
    description: "Chemical disposal table field",
    defaultResponse: [],
  },
  [INPUT_TYPES.SortableList]: {
    label: "Sortable List",
    icon: FaSort,
    description: "Sortable list field",
    defaultResponse: [],
  },
  [INPUT_TYPES.ProjectGroupPlanTable]: {
    label: "Project Group Plan Table",
    icon: FaTable,
    description: "Project group plan table field",
    defaultResponse: [],
  },
  [INPUT_TYPES.ProjectGroupHazardTable]: {
    label: "Project Group Hazard Table",
    icon: FaExclamationTriangle,
    description: "Project group hazard table field",
    defaultResponse: [],
  },
  [INPUT_TYPES.ImageFile]: {
    label: "Image File",
    icon: FaFileUpload,
    description: "Image upload field",
    defaultResponse: [],
  },
  [INPUT_TYPES.YieldTable]: {
    label: "Yield Table",
    icon: FaChartLine,
    description: "Yield calculation table field",
    defaultResponse: [],
  },
  [INPUT_TYPES.MultiYieldTable]: {
    label: "Multi Yield Table",
    icon: FaChartLine,
    description: "Multi yield calculation table field",
    defaultResponse: [],
  },
  [INPUT_TYPES.GreenMetricsTable]: {
    label: "Green Metrics Table",
    icon: FaChartBar,
    description: "Green metrics calculation table field",
    defaultResponse: [],
  },
  [INPUT_TYPES.MultiGreenMetricsTable]: {
    label: "Multi Green Metrics Table",
    icon: FaChartBar,
    description: "Multi green metrics calculation table field",
    defaultResponse: [],
  },
  [INPUT_TYPES.DateAndTime]: {
    label: "Date and Time",
    icon: FaCalendarAlt,
    description: "Date and time input field",
    defaultResponse: "",
  },
  [INPUT_TYPES.FormattedTextInput]: {
    label: "Formatted Text Input",
    icon: FaEdit,
    description: "Formatted text input field",
    defaultResponse: "",
  },
};

export const InputTypePalette = ({ onAdd }) => {
  const { data: inputTypes } = useInputTypes();
  return (
    <VStack
      px={2}
      py={4}
      align="start"
      spacing={4}
      borderWidth={1}
      borderRadius={4}
      borderColor="teal.200"
      w="full"
      maxW="220px"
    >
      <Badge label="Input Types" colorScheme="teal" />
      <Divider />
      <VStack spacing={4} align="stretch" w="full">
        {inputTypes
          .filter((inputType) => INPUT_TYPES_MAP[inputType.title])
          .map((inputType) => {
            return (
              <InputTypeItem
                key={inputType.id}
                inputType={inputType}
                onAdd={onAdd}
              />
            );
          })}
      </VStack>
    </VStack>
  );
};

const InputTypeItem = ({ inputType, onAdd }) => {
  const dragRef = useRef(null);
  const [isGrabbed, setIsGrabbed] = useState(false);

  useEffect(() => {
    if (!dragRef.current) return;

    const cleanup = draggable({
      element: dragRef.current,
      getInitialData: () => ({ inputType, type: DRAG_TYPES.INPUT_TYPE }),

      onDragEnd: () => {
        setIsGrabbed(false);
      },
    });

    return cleanup;
  }, [inputType]);

  const handleDoubleClick = () => {
    onAdd(inputType);
  };

  return (
    <HStack
      w="100%"
      ref={dragRef}
      borderRadius="xl"
      spacing={2}
      borderWidth={1}
      borderColor="gray.200"
      py={0.2}
      px={2}
      cursor={isGrabbed ? "grabbing" : "grab"}
      onMouseDown={() => setIsGrabbed(true)}
      onMouseUp={() => setIsGrabbed(false)}
      onMouseLeave={() => setIsGrabbed(false)}
      onDoubleClick={handleDoubleClick}
      _hover={{
        bg: "teal.50",
        borderColor: "teal.300",
        cursor: isGrabbed ? "grabbing" : "grab",
      }}
    >
      <Icon as={INPUT_TYPES_MAP[inputType.title].icon} fontSize="sm" />
      <Text fontSize="xs" fontWeight="normal" flex={1}>
        {INPUT_TYPES_MAP[inputType.title].label}
      </Text>
    </HStack>
  );
};
