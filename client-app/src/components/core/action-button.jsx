import {
  Box,
  Button,
  IconButton,
  Menu,
  MenuButton,
  MenuItem,
  MenuList,
  Portal,
  Text,
} from "@chakra-ui/react";
import { useMemo } from "react";
import { FaBattleNet, FaChevronDown, FaEllipsisH } from "react-icons/fa";

export const ActionButton = ({
  label,
  LeftIcon = FaBattleNet,
  actions,
  ...p
}) => {
  const eligibleActions = useMemo(
    () => Object.values(actions).filter((action) => action.isEligible()),
    [actions],
  );

  return (
    <Box>
      <Menu>
        <MenuButton
          as={label ? Button : IconButton}
          leftIcon={label && <LeftIcon />}
          rightIcon={label && <FaChevronDown />}
          aria-label={label ? label : "Actions"}
          icon={<FaEllipsisH />}
          variant={label ? "outline" : "ghost"}
          colorScheme="gray"
          size="xs"
          {...p}
        >
          <Text fontWeight="medium">{label}</Text>
        </MenuButton>
        <Portal>
          <MenuList>
            {eligibleActions.length === 0 ? (
              <MenuItem isDisabled>No actions available.</MenuItem>
            ) : (
              eligibleActions.map((x, index) => (
                <MenuItem
                  icon={x.icon}
                  onClick={x.onClick}
                  key={index}
                  isDisabled={x.isDisabled}
                >
                  {x.label}
                </MenuItem>
              ))
            )}
          </MenuList>
        </Portal>
      </Menu>
    </Box>
  );
};
