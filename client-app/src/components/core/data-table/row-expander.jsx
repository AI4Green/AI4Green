import { IconButton } from "@chakra-ui/react";
import { FaChevronDown, FaChevronRight } from "react-icons/fa";

export const DataTableRowExpander = ({ row }) =>
  row.getCanExpand() ? (
    <IconButton
      size="xs"
      icon={row.getIsExpanded() ? <FaChevronRight /> : <FaChevronDown />}
      variant="ghost"
      onClick={() => row.toggleExpanded()}
      paddingLeft={row.depth * 3}
    />
  ) : null;
