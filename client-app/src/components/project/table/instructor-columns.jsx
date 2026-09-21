import { Avatar, HStack, Text } from "@chakra-ui/react";
import { ActionButton } from "components/core/action-button";
import { Badge } from "components/core/Badge";
import { DataTableColumnHeader } from "components/core/data-table";
import { PROJECTMANAGEMENT_PERMISSIONS } from "constants";
import { useUser } from "contexts";
import { FaTrash } from "react-icons/fa";
import { useSearchParams } from "react-router-dom";

export const instructorColumns = [
  {
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Name" />
    ),
    accessorKey: "name",
    cell: ({ cell }) => (
      <HStack>
        <Avatar size="sm" name={cell.getValue() || "⚠️"} />
        <Text as="b">{cell.getValue()}</Text>
      </HStack>
    ),
    enableHiding: false,
  },
  {
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Email" />
    ),
    accessorKey: "email",
    cell: ({ cell, row }) => (
      <Text
        as="b"
        color={row.original.emailConfirmed ? "green.500" : "orange.500"}
      >
        {cell.getValue()}
      </Text>
    ),
  },
  {
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Roles" />
    ),
    accessorKey: "roles",
    cell: ({ cell }) => {
      const roles = cell.getValue();
      return (
        cell.getValue()?.length > 0 && (
          <HStack>
            {roles.map((role) => (
              <Badge
                key={role}
                colorScheme="blue"
                label={role}
                variant="outline"
                fontSize="xxs"
              />
            ))}
          </HStack>
        )
      );
    },
  },
  {
    id: "actions",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Actions" />
    ),
    cell: ({ row }) => {
      return <Actions instructor={row.original} />;
    },
    maxSize: 5,
  },
];

const Actions = ({ instructor }) => {
  const { user } = useUser();
  const [, setSearchParams] = useSearchParams();

  const actions = {
    remove: {
      isEligible: () =>
        user?.permissions?.includes(
          PROJECTMANAGEMENT_PERMISSIONS.InviteInstructors,
        ),
      icon: <FaTrash />,
      label: "Remove",
      onClick: () =>
        setSearchParams({
          action: "remove-instructor",
          instructorId: instructor.id,
        }),
      colorScheme: "red",
    },
  };
  return <ActionButton actions={actions} size="xs" />;
};
