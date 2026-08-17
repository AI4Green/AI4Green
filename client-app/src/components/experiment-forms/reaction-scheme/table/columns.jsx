import { Text } from "@chakra-ui/react";
import {
  DataTableColumnHeader,
  TableCellCheckBox,
  TableCellDeleteRowButton,
  TableCellDropdown,
  TableCellNumberInput,
  TableCellNumberInputWithUnit,
} from "components/core/data-table";

import { HazardsInput } from "./reaction-table";

/**
 *
 * @param {*} config - contains the following properties:
 * - isDisabled: boolean (whether the table is disabled or not)
 * - massColumnHeaderDropdown: object (contains the following properties)
 *    - ColumnUnitHeaderDropdown: dropdown component
 *    - props: props to pass to the dropdown
 * @returns - an array of objects to be used as columns in the DataTable
 */

export const reactionTableColumns = ({ isDisabled }) => [
  {
    accessorKey: "substanceType",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Type" />
    ),
    cell: ({ getValue }) => <Text fontSize="sm">{getValue()}</Text>,
  },
  {
    accessorKey: "substancesUsed",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Substances Used" />
    ),
  },
  {
    accessorKey: "limiting",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Limiting" />
    ),
    cell: ({ getValue, row, column, table }) => (
      <TableCellCheckBox
        {...{ getValue, row, column, table }}
        isDisabled={isDisabled}
      />
    ),
  },
  {
    accessorKey: "mass",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Mass (Vol)" />
    ),
    cell: ({ getValue, row, column, table }) => (
      <TableCellNumberInputWithUnit
        {...{ getValue, row, column, table }}
        isDisabled={isDisabled}
        placeholder="Mass (Vol)"
        options={[
          { value: "cm³", label: "cm³" },
          { value: "g", label: "g" },
        ]}
      />
    ),
  },
  {
    accessorKey: "gls",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="g/l/s (Physical form)" />
    ),
    cell: ({ getValue, row, column, table }) => (
      <TableCellDropdown
        {...{ getValue, row, column, table }}
        options={[
          { value: "g", label: "Gas" },
          { value: "l", label: "Liquid" },
          { value: "s", label: "Solid" },
        ]}
        isDisabled={isDisabled}
      />
    ),
  },
  {
    accessorKey: "molWeight",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Mol.Wt" />
    ),
  },
  {
    accessorKey: "amount",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Amount (mmol)" />
    ),
    cell: ({ getValue, row, column, table }) => (
      <TableCellNumberInput
        {...{ getValue, row, column, table }}
        isDisabled={isDisabled}
        placeholder="Amount (mmol)"
      />
    ),
  },
  {
    accessorKey: "density",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Density" />
    ),
  },
  {
    accessorKey: "hazardsInput",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Hazards" w={48} />
    ),
    cell: ({ getValue, row, column, table }) => (
      <HazardsInput
        {...{ getValue, row, column, table }}
        isDisabled={isDisabled}
      />
    ),
  },
  ...(!isDisabled
    ? [
        {
          id: "actions",
          header: ({ column }) => (
            <DataTableColumnHeader column={column} title="Actions" />
          ),
          cell: ({ row, table }) =>
            row.original?.manualEntry ? (
              <TableCellDeleteRowButton {...{ row, table }} /> // only show delete button if the row is manually entered
            ) : null,
        },
      ]
    : []),
];
