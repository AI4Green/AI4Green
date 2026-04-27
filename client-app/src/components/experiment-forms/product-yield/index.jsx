import { Button, HStack, Text, VStack } from "@chakra-ui/react";
import { DataTable } from "components/core/data-table";
import { useField } from "formik";
import { useEffect, useMemo, useState } from "react";
import { FaPlus } from "react-icons/fa";

import { productYieldTableColumn } from "./columns";

export const ProductYieldTable = ({ name, label, isDisabled }) => {
  const [field, , helpers] = useField(name);
  const [tableData, setTableData] = useState(field.value || []);

  useEffect(() => {
    if (tableData !== field.value) {
      helpers.setValue(tableData);
    }
  }, [tableData, field.value, helpers]);

  const columns = useMemo(
    () => productYieldTableColumn(isDisabled),
    [isDisabled],
  );

  return (
    <VStack align="flex-start" w="full">
      <DataTable
        data={tableData}
        setTableData={setTableData}
        columns={columns}
        FooterCellAddRow={
          !isDisabled && (
            <AddRowButton
              handleAddRow={() =>
                handleAddRow(columns, setTableData, tableData)
              }
            />
          )
        }
      >
        <HStack flex={1}>
          <Text size="sm" as="b">
            {label}
          </Text>
        </HStack>
      </DataTable>
    </VStack>
  );
};

const AddRowButton = ({ handleAddRow }) => {
  return (
    <Button
      colorScheme="blue"
      size="xs"
      leftIcon={<FaPlus />}
      onClick={handleAddRow}
    >
      Add new
    </Button>
  );
};

const handleAddRow = (columns, setTableData, tableData) => {
  const defaultValues = {
    serialNumber: tableData.length + 1,
    product: "",
    expectedMass: { unit: "", value: 0 },
    actualMass: { unit: "", value: 0 },
    moles: 0,
    yield: 0,
  };

  const newRow = columns
    .filter((column) => "accessorKey" in column)
    .reduce((acc, column) => {
      return {
        ...acc,
        [column.accessorKey]: defaultValues[column.accessorKey] ?? "",
      };
    }, {});

  setTableData((old) => [...old, newRow]);
};
