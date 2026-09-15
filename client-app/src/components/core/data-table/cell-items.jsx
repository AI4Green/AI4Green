import {
  Badge,
  Box,
  Checkbox,
  IconButton,
  Input,
  NumberDecrementStepper,
  NumberIncrementStepper,
  NumberInput,
  NumberInputField,
  NumberInputStepper,
  Select,
  Text,
  Textarea,
  VStack,
} from "@chakra-ui/react";
import { countWords } from "helpers";
import { useCallback, useEffect, useState } from "react";
import { FaRegTimesCircle } from "react-icons/fa";
import { MdArrowDropDown } from "react-icons/md";

const COMPONENT_TYPES = {
  TextInput: "TextInput",
  TextAreaInput: "TextAreaInput",
  DateInput: "DateInput",
  Dropdown: "Dropdown",
  NumberInput: "NumberInput",
  CheckBox: "CheckBox",
};

const getInitialValue = (componentType, getValue) => {
  const { NumberInput, CheckBox } = COMPONENT_TYPES;
  switch (componentType) {
    case NumberInput:
      return getValue() || 0;
    case CheckBox:
      return getValue() || false;
    default:
      return getValue() || "";
  }
};

const getOnChangeHandler = (componentType, setValue, value) => {
  const { NumberInput, CheckBox } = COMPONENT_TYPES;
  switch (componentType) {
    case NumberInput:
      return (e) => setValue(parseFloat(e));
    case CheckBox:
      return () => setValue(!value);
    default:
      return (e) => setValue(e.target.value);
  }
};

const withTableCell = (Component, componentType) => {
  const WrappedComponent = ({ getValue, row, column, table, ...props }) => {
    const initialValue = getInitialValue(componentType, getValue);
    const [value, setValue] = useState(initialValue);

    const onBlur = useCallback(() => {
      table.options.meta?.updateData(row.index, column.id, value);
    }, [table.options.meta, row.index, column.id, value]);

    const onChange = getOnChangeHandler(componentType, setValue, value);

    useEffect(() => {
      setValue(initialValue);
    }, [initialValue]);

    const componentProps =
      componentType === COMPONENT_TYPES.CheckBox
        ? { isChecked: value, onChange, onBlur, ...props }
        : { value, onChange: onChange, onBlur, ...props };

    return <Component {...componentProps} />;
  };

  WrappedComponent.displayName = `withTableCell(${
    Component.displayName || Component.name || "Component"
  })`;
  return WrappedComponent;
};

export const TableCellTextInput = withTableCell(
  (props) => <Input size="sm" {...props} />,
  COMPONENT_TYPES.TextInput,
);

export const TableCellTextAreaInput = withTableCell(
  (props) => (
    <VStack align="flex-e">
      <Textarea size="sm" rows="6" {...props} _disabled={{ opacity: 0.7 }} />
      <WordCountBadge {...props} />
    </VStack>
  ),
  COMPONENT_TYPES.TextAreaInput,
);

export const TableCellDateInput = withTableCell(
  (props) => (
    <Input size="sm" type="date" {...props} _disabled={{ opacity: 0.7 }} />
  ),
  COMPONENT_TYPES.DateInput,
);

export const TableCellDropdown = withTableCell(
  ({ options, ...props }) => (
    <Box maxW="130px">
      <Select size="sm" placeholder="Select option" {...props}>
        {options.map((option, index) => (
          <option key={index} value={option.value}>
            {option.label}
          </option>
        ))}
      </Select>
    </Box>
  ),
  COMPONENT_TYPES.Dropdown,
);

export const TableCellNumberInput = withTableCell(
  (props) => (
    <Box maxW="100px">
      <NumberInput size="sm" step={0.2} {...props} _disabled={{ opacity: 0.7 }}>
        <NumberInputField borderRadius={4} />
        <NumberInputStepper>
          <NumberIncrementStepper />
          <NumberDecrementStepper />
        </NumberInputStepper>
      </NumberInput>
    </Box>
  ),
  COMPONENT_TYPES.NumberInput,
);

export const TableCellCheckBox = withTableCell(
  (props) => <Checkbox {...props} />,
  COMPONENT_TYPES.CheckBox,
);

export const TableCellDeleteRowButton = ({ row, table }) => (
  <IconButton
    onClick={() => table.options.meta?.removeRow(row.index)}
    aria-label="Delete Row"
    icon={<FaRegTimesCircle />}
    variant="ghost"
    fontSize="md"
    colorScheme="red"
  />
);

export const WordCountBadge = ({ value, limit }) => {
  const count = countWords(value);
  return limit ? (
    <Badge colorScheme={count > limit ? "red" : undefined}>
      Word Count: {count} / {limit}
    </Badge>
  ) : null;
};

export const TableCellNumberInputWithUnit = ({
  getValue,
  row,
  column,
  table,
  placeholder,
  isDisabled,
  options = [],
}) => {
  const initialValue = getValue();
  const [unit, setUnit] = useState(initialValue?.unit);
  const [value, setValue] = useState(initialValue?.value);
  const [error, setError] = useState("");

  const onBlur = useCallback(() => {
    if (value && value !== 0 && !unit) {
      setError("Please select a unit");
      return;
    }

    setError("");

    table.options.meta?.updateData(row.index, column.id, {
      unit,
      value,
    });
  }, [table, row.index, column.id, unit, value]);

  useEffect(() => {
    setUnit(initialValue?.unit);
    setValue(initialValue?.value);
  }, [initialValue]);

  useEffect(() => {
    if (!unit) {
      setValue(0);
    }
  }, [unit]);

  return (
    <VStack align="center" maxW={24}>
      <NumberInput
        size="sm"
        step={0.2}
        onChange={(v) => setValue(v)}
        onBlur={onBlur}
        value={value}
        isDisabled={isDisabled}
        placeholder={placeholder}
        _disabled={{ opacity: 0.7 }}
      >
        <NumberInputField borderRadius={4} />
        <NumberInputStepper>
          <NumberIncrementStepper />
          <NumberDecrementStepper />
        </NumberInputStepper>
      </NumberInput>

      <Select
        icon={<MdArrowDropDown />}
        size="xs"
        placeholder="Unit"
        onChange={(e) => setUnit(e.target.value)}
        onBlur={onBlur}
        isDisabled={isDisabled}
        borderRadius={4}
        value={unit}
      >
        {options.map((option) => (
          <option key={option.value} value={option.value}>
            {option.label}
          </option>
        ))}
      </Select>
      {error && <Text color="red">{error}</Text>}
    </VStack>
  );
};
