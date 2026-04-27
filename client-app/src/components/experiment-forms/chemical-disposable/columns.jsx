import {
  DataTableColumnHeader,
  TableCellCheckBox,
  TableCellDeleteRowButton,
  TableCellTextInput,
} from "components/core/data-table";

import { TableCellOther } from ".";

export const chemicalDisposableTableColumns = ({ isDisabled }) => [
  {
    accessorKey: "serialNumber",
    header: () => <DataTableColumnHeader sorting={false} title="No." />,
    enableHiding: false,
  },
  {
    accessorKey: "chemical",
    header: () => <DataTableColumnHeader sorting={false} title="Chemical" />,
    cell: ({ getValue, row, column, table }) => (
      <TableCellTextInput
        {...{ getValue, row, column, table }}
        isDisabled={isDisabled}
        placeholder="Chemical"
      />
    ),
  },
  {
    accessorKey: "halogenatedSolvent",
    header: () => (
      <DataTableColumnHeader sorting={false} title="Halogenated Solvent" />
    ),
    cell: ({ getValue, row, column, table }) => (
      <TableCellCheckBox
        {...{ getValue, row, column, table }}
        isDisabled={isDisabled}
      />
    ),
  },
  {
    accessorKey: "nonHalogenatedSolvent",
    header: () => (
      <DataTableColumnHeader sorting={false} title="Non-Halogenated Solvent" />
    ),
    cell: ({ getValue, row, column, table }) => (
      <TableCellCheckBox
        {...{ getValue, row, column, table }}
        isDisabled={isDisabled}
      />
    ),
  },
  {
    accessorKey: "wasteStore",
    header: () => <DataTableColumnHeader sorting={false} title="Waste Store" />,
    cell: ({ getValue, row, column, table }) => (
      <TableCellCheckBox
        {...{ getValue, row, column, table }}
        isDisabled={isDisabled}
      />
    ),
  },
  {
    accessorKey: "aqueousNonHazardous",
    header: () => (
      <DataTableColumnHeader sorting={false} title="Aqueous (Non Hazardous)" />
    ),
    cell: ({ getValue, row, column, table }) => (
      <TableCellCheckBox
        {...{ getValue, row, column, table }}
        isDisabled={isDisabled}
      />
    ),
  },
  {
    accessorKey: "others",
    header: () => (
      <DataTableColumnHeader sorting={false} title="Other (specify)" />
    ),
    cell: ({ getValue, row, column, table }) => (
      <TableCellOther
        {...{ getValue, row, column, table }}
        isDisabled={isDisabled}
      />
    ),
  },
  ...(!isDisabled
    ? [
        {
          id: "actions",
          header: () => (
            <DataTableColumnHeader sorting={false} title="Actions" />
          ),
          cell: TableCellDeleteRowButton,
        },
      ]
    : []),
];
