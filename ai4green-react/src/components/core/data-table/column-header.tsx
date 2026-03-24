import {
  Box,
  Button,
  Icon,
  Menu,
  MenuButton,
  MenuDivider,
  MenuItem,
  MenuList,
  Text,
} from "@chakra-ui/react";
import { FaChevronDown, FaChevronUp, FaEyeSlash, FaSort } from "react-icons/fa";

export const DataTableColumnHeader = ({ sorting, column, title, w = 0 }) => {
  if (!sorting && !column?.getCanSort()) {
    return (
      <Text
        py={2}
        fontSize="xs"
        fontWeight="medium"
        noOfLines={3}
        color="gray.800"
        letterSpacing="tight"
      >
        {title}
      </Text>
    );
  }

  return (
    <Box align="start" w={w}>
      <Menu>
        <MenuButton
          as={Button}
          variant="ghost"
          size="xs"
          px={0}
          rightIcon={
            <Icon
              as={
                column.getIsSorted() === "desc"
                  ? FaChevronDown
                  : column.getIsSorted() === "asc"
                  ? FaChevronUp
                  : FaSort
              }
              fontSize="xs"
              color={
                column.getIsSorted() === "desc"
                  ? "green"
                  : column.getIsSorted() === "asc"
                  ? "green"
                  : "gray.600"
              }
            />
          }
          _hover={{ bg: "transparent" }}
          _active={{ bg: "transparent" }}
        >
          <Text fontSize="xs" fontWeight="medium" noOfLines={3}>
            {title}
          </Text>
        </MenuButton>
        <MenuList>
          <MenuItem
            onClick={() => column.toggleSorting(false)}
            icon={<FaChevronUp />}
          >
            Asc
          </MenuItem>
          <MenuItem
            onClick={() => column.toggleSorting(true)}
            icon={<FaChevronDown />}
          >
            Desc
          </MenuItem>
          <MenuDivider />

          <MenuItem
            onClick={() => column.toggleVisibility(false)}
            icon={<FaEyeSlash />}
          >
            Hide
          </MenuItem>
        </MenuList>
      </Menu>
    </Box>
  );
};
