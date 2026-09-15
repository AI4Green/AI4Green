import {
  Box,
  Button,
  Menu,
  MenuButton,
  MenuGroup,
  MenuItem,
  MenuList,
} from "@chakra-ui/react";
import { FaCheck, FaChevronDown, FaEyeSlash } from "react-icons/fa";

export function DataTableViewOptions({ table }) {
  return (
    <Box>
      <Menu>
        <MenuButton
          as={Button}
          size="xs"
          variant="outline"
          rightIcon={<FaChevronDown />}
        >
          View
        </MenuButton>
        <MenuList align="end" minWidth="150px">
          <MenuGroup title="Table view options">
            {table
              .getAllColumns()
              .filter(
                (column) =>
                  typeof column.accessorFn !== "undefined" &&
                  column.getCanHide(),
              )
              .map((column) => {
                return (
                  <MenuItem
                    icon={column.getIsVisible() ? <FaCheck /> : <FaEyeSlash />}
                    key={column.id}
                    textTransform="capitalize"
                    onClick={() =>
                      column.toggleVisibility(!column.getIsVisible())
                    }
                  >
                    {column.id}
                  </MenuItem>
                );
              })}
          </MenuGroup>
        </MenuList>
      </Menu>
    </Box>
  );
}
