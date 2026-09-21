import {
  DataTableColumnHeader,
  TableCellDeleteRowButton,
  TableCellTextAreaInput,
  TableCellTextInput,
} from "components/core/data-table";

export const hazardSummaryTableColumn = (isDisabled) => [
  {
    accessorKey: "serialNumber",
    header: () => <DataTableColumnHeader sorting={false} title="No." />,
    enableHiding: false,
  },
  {
    accessorKey: "chemicals",
    header: () => <DataTableColumnHeader sorting={false} title="Chemicals" />,
    cell: ({ getValue, row, column, table }) => (
      <TableCellTextInput
        {...{ getValue, row, column, table }}
        placeholder="Chemical"
        isDisabled={isDisabled}
      />
    ),
  },
  {
    accessorKey: "hazard",
    header: () => <DataTableColumnHeader sorting={false} title="Hazards" />,
    cell: ({ getValue, row, column, table }) => (
      <TableCellTextInput
        {...{ getValue, row, column, table }}
        placeholder="Hazard"
        isDisabled={isDisabled}
      />
    ),
  },
  {
    accessorKey: "precautions",
    header: () => <DataTableColumnHeader sorting={false} title="Precautions" />,
    cell: ({ getValue, row, column, table }) => (
      <TableCellTextAreaInput
        {...{ getValue, row, column, table }}
        placeholder="Precaution"
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
