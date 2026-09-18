/**
 * The component renders a table.
 * Populates table using ketcher data if there's no existing table data.
 * Data is reset if ketcher user resets the ketcher.
 */

import {
  Button,
  HStack,
  Icon,
  IconButton,
  Input,
  Popover,
  PopoverArrow,
  PopoverBody,
  PopoverContent,
  PopoverTrigger,
  Text,
  Textarea,
  useDisclosure,
  VStack,
} from "@chakra-ui/react";
import { DataTable } from "components/core/data-table";
import { useCallback, useMemo, useState } from "react";
import {
  FaCheckCircle,
  FaExclamationCircle,
  FaFlask,
  FaInfo,
  FaVial,
} from "react-icons/fa";

import { AddSubstanceModal } from "../modal/add-substance";
import { reactionTableColumns } from "./columns";
import { useReactionTable } from "./use-data";

export const REACTION_TABLE_DEFAULT_VALUES = {
  tableData: [],
  massUnit: "",
};

export const ReactionTable = ({ name, ketcherData, isDisabled }) => {
  const { tableData, setTableData } = useReactionTable(name, ketcherData);

  const tableColumns = useMemo(
    () =>
      reactionTableColumns({
        isDisabled,
      }),
    [isDisabled],
  ); // array of objects to be used for the DataTable columns

  const footerCellProps = useCallback(() => ({ setTableData }), [setTableData]);

  return (
    <DataTable
      data={tableData}
      setTableData={setTableData}
      columns={tableColumns}
      FooterCellAddRow={!isDisabled && <FooterCell {...footerCellProps()} />}
    >
      <HStack flex={1}>
        <Text size="sm" as="b">
          Please fill in the relevant fields below
        </Text>
      </HStack>
    </DataTable>
  );
};

export const FooterCell = ({ setTableData }) => {
  const Reagent = useDisclosure();
  const Solvent = useDisclosure();

  return (
    <HStack justify="flex-end" spacing={5} py={2}>
      <Button
        colorScheme="pink"
        size="sm"
        leftIcon={<FaFlask />}
        onClick={Reagent.onOpen}
      >
        Add reagent
      </Button>
      <Button
        colorScheme="teal"
        size="sm"
        leftIcon={<FaVial />}
        onClick={Solvent.onOpen}
      >
        Add solvent
      </Button>
      <AddSubstanceModal
        isModalOpen={Reagent.isOpen}
        onModalClose={Reagent.onClose}
        setTableData={setTableData}
      />
      <AddSubstanceModal
        isModalOpen={Solvent.isOpen}
        onModalClose={Solvent.onClose}
        setTableData={setTableData}
        isAddingSolvent
      />
    </HStack>
  );
};

//TODO: add validation for the hazards against the backend
const HazardsValidation = ({ input, valid }) => {
  if (!input || !valid) return null;

  return input === valid ? (
    <Icon as={FaCheckCircle} color="green.500" fontSize="xs" mx={1} />
  ) : (
    <HoverContent
      ariaLabelIcon="warning"
      ariaLabelContent="hazard-warning"
      Icon={FaExclamationCircle}
      color="orange.500"
      content="Incorrect hazard code."
    />
  );
};

export const HazardsInput = ({ getValue, row, column, table, isDisabled }) => {
  const initialValue = getValue();
  const [hazardCodes, setHazardCodes] = useState(
    initialValue?.hazardCodes || "",
  );
  const [hazardDescription, setHazardDescription] = useState(
    initialValue?.hazardDescription || "",
  );

  const onBlur = () => {
    table.options.meta?.updateData(row.index, column.id, {
      hazardCodes,
      hazardDescription,
    });
  };

  return (
    <VStack align="stretch" spacing={3}>
      <HStack w="full" spacing={1}>
        <HazardsValidation input={hazardCodes} valid={row.original.hazards} />
        <Input
          size="sm"
          value={hazardCodes}
          onChange={(e) => setHazardCodes(e.target.value)}
          onBlur={onBlur}
          placeholder="Hazard Codes"
          isDisabled={isDisabled}
          flex={1}
        />
        <HoverContent
          ariaLabelIcon="information"
          ariaLabelContent="hazard code input example"
          Icon={FaInfo}
          color="gray.500"
          content="For example, the hazard codes for Methane are H220, H280, and H281. Enter them as H220-H280-H281."
        />
      </HStack>
      <Textarea
        size="xs"
        value={hazardDescription}
        onChange={(e) => setHazardDescription(e.target.value)}
        onBlur={onBlur}
        placeholder="Hazard Codes Description"
        isDisabled={isDisabled}
        rows={2}
      />
    </VStack>
  );
};

const HoverContent = ({
  ariaLabelIcon,
  ariaLabelContent,
  Icon,
  color,
  contentColor = "white",
  content,
}) => (
  <Popover>
    <PopoverTrigger>
      <IconButton
        aria-label={ariaLabelIcon}
        icon={<Icon />}
        size="xs"
        variant="ghost"
        color={color}
      />
    </PopoverTrigger>
    <PopoverContent
      color={contentColor}
      bg={color}
      aria-label={ariaLabelContent}
      w="auto"
      maxW="xs"
    >
      <PopoverArrow />
      <PopoverBody p={2} border="0">
        <Text fontSize="xs" whiteSpace="normal">
          {content}
        </Text>
      </PopoverBody>
    </PopoverContent>
  </Popover>
);
