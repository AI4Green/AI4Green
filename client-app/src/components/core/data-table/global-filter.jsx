import { HStack, Input, InputGroup, InputLeftElement } from "@chakra-ui/react";
import { FaSearch } from "react-icons/fa";

export const DataTableGlobalFilter = ({
  setSearchValue,
  searchValue,
  placeholder,
}) => (
  <HStack>
    <InputGroup>
      <InputLeftElement pointerEvents="none" height="100%">
        <FaSearch color="gray" />
      </InputLeftElement>
      <Input
        size="sm"
        placeholder={placeholder}
        _placeholder={{ opacity: 1 }}
        onChange={(e) => setSearchValue(e.target.value)}
        value={searchValue || ""}
      />
    </InputGroup>
  </HStack>
);
