import {
  Box,
  HStack,
  Table,
  TableContainer,
  Tbody,
  Td,
  Tfoot,
  Th,
  Thead,
  Tr,
} from "@chakra-ui/react";
import {
  flexRender,
  getCoreRowModel,
  getExpandedRowModel,
  getFilteredRowModel,
  getSortedRowModel,
  useReactTable,
} from "@tanstack/react-table";
import { useState } from "react";
import { useNavigate } from "react-router-dom";

import { DataTableViewOptions } from "./view-options";

/**
 * DataTable component for rendering table view.
 * Props:
 * - columns: Array of column definitions. If row is expandable, add 'expander' as column id.
 * - data: Contains the rows to be displayed in the table. Each object in the array represents a row.
 *    For expandable rows, include a 'subRows' property with an array of sub-row data.
 *    For clickable rows, include a 'targetPath' property specifying the navigation path.
 * - setTableData: Function to update the data.
 * - globalFilter: Global filter value.
 * - FooterCellAddRow: Footer cell for adding new row.
 * - children: Additional components to be rendered.
 */
export function DataTable({
  columns,
  data,
  setTableData,
  globalFilter,
  FooterCellAddRow, // TODO: In future, we can add make this flexible.
  children,
}) {
  const [sorting, setSorting] = useState([]);
  const [expanded, setExpanded] = useState({});
  const navigate = useNavigate();

  const handleRowClick = (path) => path && navigate(path);

  const handleCellClick = (e, cell, row) => {
    if (!row.original.targetPath) return;

    const isEmptyExpander =
      cell.column.id === "expander" && !row.getCanExpand();
    const isContentNotClickable =
      ["string", "number"].includes(typeof cell.getValue()) || isEmptyExpander;

    if (!isContentNotClickable) e.stopPropagation(); // Prevent row click
  };

  const table = useReactTable({
    data,
    columns,
    state: {
      sorting,
      expanded,
      globalFilter,
    },
    meta: {
      updateData: (rowIndex, columnId, value) =>
        setTableData((prev) =>
          prev.map((row, index) =>
            index === rowIndex
              ? {
                  ...prev[rowIndex],
                  [columnId]: value,
                }
              : row,
          ),
        ),
      removeRow: (rowIndex) => {
        setTableData((prev) =>
          prev
            .filter((_row, index) => index !== rowIndex)
            .map((row, index) => ({ ...row, serialNumber: index + 1 })),
        );
      },
    },
    getSubRows: (row) => row.subRows,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
    getExpandedRowModel: getExpandedRowModel(),
    onSortingChange: setSorting,
    onExpandedChange: setExpanded,
    filterFromLeafRows: true,
  });

  return (
    <Box w="full">
      <HStack justify="flex-end" py={2} mb={2} spacing={4}>
        {children}
        <DataTableViewOptions table={table} />
      </HStack>
      <TableContainer borderRadius={4} borderWidth={1}>
        <Table variant="simple" colorScheme="gray" size="sm">
          <Thead bgColor="gray.50">
            {table.getHeaderGroups().map((headerGroup) => (
              <Tr key={headerGroup.id}>
                {headerGroup.headers.map((header) => (
                  <Th
                    key={header.id}
                    textTransform="none"
                    whiteSpace="normal"
                    w={header.column.columnDef.meta?.width || "auto"}
                    py={2}
                    borderBottom="1px"
                    borderBottomColor="brand.200"
                  >
                    {header.isPlaceholder
                      ? null
                      : flexRender(
                          header.column.columnDef.header,
                          header.getContext(),
                        )}
                  </Th>
                ))}
              </Tr>
            ))}
          </Thead>
          <Tbody>
            {table.getRowModel().rows?.length ? (
              table.getRowModel().rows.map((row) => (
                <Tr
                  tabIndex={0}
                  key={row.id}
                  data-state={row.getIsSelected() && "selected"}
                  bgColor={row.getIsExpanded() ? "blackAlpha.100" : ""}
                  _hover={{
                    bg: "blackAlpha.50",
                    cursor: row.original.targetPath ? "pointer" : "default",
                  }}
                  onClick={() => handleRowClick(row.original.targetPath)}
                  onKeyDown={(e) => {
                    if (e.key === "Enter") {
                      handleRowClick(row.original.targetPath);
                    }
                  }}
                  transition="all 0.3s ease-in-out"
                >
                  {row.getVisibleCells().map((cell) => (
                    <Td
                      key={cell.id}
                      onClick={(e) => handleCellClick(e, cell, row)}
                      py={4}
                      verticalAlign={
                        cell.column.columnDef.meta?.verticalAlign || "middle"
                      }
                      maxW={cell.column.columnDef.meta?.width || "auto"}
                      overflow="hidden"
                      whiteSpace="normal"
                      wordBreak="break-words"
                    >
                      {flexRender(
                        cell.column.columnDef.cell,
                        cell.getContext(),
                      )}
                    </Td>
                  ))}
                </Tr>
              ))
            ) : (
              <Tr>
                <Td colSpan={columns.length} textAlign="center">
                  No data available
                </Td>
              </Tr>
            )}
          </Tbody>
          {FooterCellAddRow && (
            <Tfoot>
              <Tr>
                <Th colSpan={columns.length} textAlign="right">
                  {FooterCellAddRow}
                </Th>
              </Tr>
            </Tfoot>
          )}
        </Table>
      </TableContainer>
    </Box>
  );
}
