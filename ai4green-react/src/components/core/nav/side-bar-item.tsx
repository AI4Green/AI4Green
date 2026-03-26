import {
  Box,
  HStack,
  Icon,
  LinkBox,
  LinkOverlay,
  Text,
} from "@chakra-ui/react";
import { Link } from "react-router-dom";

export const SidebarButton = ({ item, onClose }) => {
  if (!item.path) {
    return <SidebarTitle item={item} />;
  }

  return <SidebarLinkButton item={item} onClose={onClose} />;
};

const SidebarLinkButton = ({ item, onClose }) => (
  <LinkBox
    w="full"
    _hover={{
      color: "brand.500",
      transition: "color 0.2s ease-in-out",
    }}
  >
    <LinkOverlay as={Link} to={item.path} onClick={onClose}>
      <HStack>
        <Icon as={item.icon} boxSize={4} marginRight={2} />
        <Text fontSize={{ base: "sm", md: "md" }} fontWeight="light" py={1}>
          {item.label}
        </Text>
      </HStack>
    </LinkOverlay>
  </LinkBox>
);

const SidebarTitle = ({ item }) => {
  return (
    <Box display="flex" alignItems="center">
      <Icon as={item.icon} boxSize={5} marginRight={2} />
      <Text fontSize="lg" fontWeight="light" letterSpacing="tight" py={1}>
        {item.label}
      </Text>
    </Box>
  );
};

export default SidebarButton;
