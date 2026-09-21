import { Button, FormLabel, HStack, VStack } from "@chakra-ui/react";
import { DataTable } from "components/core/data-table";
import { useFormikContext } from "formik";
import { useEffect, useMemo, useState } from "react";
import { FaPlus } from "react-icons/fa";

import { hazardSummaryTableColumn } from "./columns";

export const HazardSummaryTable = ({ name, label, isDisabled }) => {
  const { values, setFieldValue } = useFormikContext();
  const [tableData, setTableData] = useState(values[name] || []);

  useEffect(() => {
    if (tableData !== values[name]) {
      setFieldValue(name, tableData);
    }
  }, [tableData, name, setFieldValue, values]);

  const columns = useMemo(
    () => hazardSummaryTableColumn(isDisabled),
    [isDisabled],
  );

  const handleAddRow = () => {
    const newRow = columns
      .filter((column) => "accessorKey" in column) // filter out columns without property 'accessorKey'
      .reduce((acc, column) => {
        return {
          ...acc,
          [column.accessorKey]:
            column.accessorKey === "serialNumber" ? tableData.length + 1 : "",
        };
      }, {});

    setTableData((old) => [...old, newRow]);
  };

  return (
    <VStack align="flex-start" w="full">
      <DataTable
        data={tableData}
        setTableData={setTableData}
        columns={columns}
        FooterCellAddRow={
          !isDisabled && <FooterCell handleAddRow={handleAddRow} />
        }
      >
        <HStack flex={1}>
          <FormLabel>{label}</FormLabel>
        </HStack>
      </DataTable>
    </VStack>
  );
};

export const FooterCell = ({ handleAddRow }) => {
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
