import { Badge as ChakraBadge, HStack, Icon, Text } from "@chakra-ui/react";

export const Badge = ({ label, leftIcon, rightIcon, colorScheme, ...p }) => {
  return (
    <ChakraBadge
      colorScheme={colorScheme}
      px={4}
      py={0.5}
      borderRadius="xl"
      fontSize="xs"
      fontWeight="medium"
      textTransform="capitalize"
      {...p}
    >
      <HStack spacing={2}>
        {leftIcon && <Icon as={leftIcon} />}
        <Text>{label}</Text>
        {rightIcon && <Icon as={rightIcon} />}
      </HStack>
    </ChakraBadge>
  );
};
