import { Flex, Text } from "@chakra-ui/react";

export default function Home() {
  return (
    <Flex bg="blue.500" p={4} color="white" align="center">
      <Text fontSize="xl" fontWeight="bold">
        My App
      </Text>
    </Flex>
  );
}
