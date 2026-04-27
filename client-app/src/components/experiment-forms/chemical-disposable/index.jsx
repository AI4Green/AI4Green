import { Button, Checkbox, FormLabel, HStack, Input } from "@chakra-ui/react";
import { DataTable } from "components/core/data-table";
import { useFormikContext } from "formik";
import { useEffect, useMemo, useState } from "react";
import { FaPlus } from "react-icons/fa";

import { chemicalDisposableTableColumns } from "./columns";

export const ChemicalDisposableTable = ({ name, label, isDisabled }) => {
  const { values, setFieldValue } = useFormikContext();
  const [tableData, setTableData] = useState(values[name] || []);

  useEffect(() => {
    if (tableData !== values[name]) {
      setFieldValue(name, tableData);
    }
  }, [tableData, name, setFieldValue, values]);

  return (
    <CDTable
      tableData={tableData}
      setTableData={setTableData}
      tableLabel={label}
      isDisabled={isDisabled}
    />
  );
};

const CDTable = ({ tableData, setTableData, tableLabel, isDisabled }) => {
  const columns = useMemo(
    () => chemicalDisposableTableColumns({ isDisabled }),
    [isDisabled],
  ); // array of objects to be used for the DataTable columns

  const handleAddRow = () => {
    const newRow = columns
      .filter((column) => "accessorKey" in column) // filter out columns without property 'accessorKey'
      .reduce((acc, column) => {
        return {
          ...acc,
          [column.accessorKey]:
            column.accessorKey === "serialNumber"
              ? tableData.length + 1
              : column.accessorKey === "chemical"
              ? ""
              : false,
        };
      }, {});

    setTableData((old) => [...old, newRow]);
  };

  return (
    <DataTable
      data={tableData}
      setTableData={setTableData}
      columns={columns}
      FooterCellAddRow={
        !isDisabled && <FooterCell handleAddRow={handleAddRow} />
      }
    >
      <HStack flex={1}>
        <FormLabel>{tableLabel}</FormLabel>
      </HStack>
    </DataTable>
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

export const TableCellOther = ({
  getValue,
  row,
  column,
  table,
  isDisabled,
}) => {
  const initialValue = getValue();
  const [status, setStatus] = useState(initialValue?.status);
  const [value, setValue] = useState(initialValue?.value);

  const onBlur = () => {
    table.options.meta?.updateData(row.index, column.id, {
      status,
      value,
    });
  };

  useEffect(() => {
    setStatus(initialValue?.status);
    setValue(initialValue?.value);
  }, [initialValue]);

  useEffect(() => {
    if (!status) {
      setValue("");
    }
  }, [status]);

  return (
    <HStack align="center">
      <Checkbox
        isChecked={status}
        onChange={() => setStatus(!status)}
        onBlur={onBlur}
        isDisabled={isDisabled}
      />

      {status && (
        <Input
          size="sm"
          value={value}
          onChange={(e) => setValue(e.target.value)}
          onBlur={onBlur}
          placeholder="Other specify"
          isDisabled={isDisabled}
        />
      )}
    </HStack>
  );
};
