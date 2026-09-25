import { Avatar, Heading, HStack, Text, VStack } from "@chakra-ui/react";
import { Badge } from "components/core/Badge";
import { TITLE_ICON_COMPONENTS } from "constants";
import { useUser } from "contexts";

export const SectionHeader = ({ header, project, owner, action }) => (
  <HStack w="full" justify="space-between" borderBottomWidth={1} py={4}>
    <VStack spacing={4} align="start">
      <HStack spacing={6}>
        {header.title && (
          <Heading
            as="h2"
            fontSize={{ base: "sm", lg: "md" }}
            fontWeight="normal"
            color="gray.700"
          >
            {header.title}
          </Heading>
        )}
      </HStack>
    </VStack>

    <VStack align="end">{action}</VStack>
  </HStack>
);

export const ExperimentHeading = ({ projectName, projectGroupName, owner }) => {
  const { user } = useUser();
  const isAuthor = owner?.id === user.userId;
  return (
    <HStack align="center" spacing={4}>
      {owner && !isAuthor && (
        <HStack>
          <Avatar name={owner.name} size="xs" />
          <Text
            fontSize={{ base: "xs", md: "sm" }}
            fontWeight="semibold"
            color="gray.700"
          >
            {owner.name}
          </Text>
        </HStack>
      )}
      <Badge
        label={projectName}
        leftIcon={TITLE_ICON_COMPONENTS.Project}
        colorScheme="brand"
        variant="outline"
        fontSize="xxs"
      />
      {projectGroupName && (
        <Badge
          label={projectGroupName}
          leftIcon={TITLE_ICON_COMPONENTS.ProjectGroup}
          colorScheme="blue"
          variant="outline"
          fontSize="xxs"
        />
      )}
    </HStack>
  );
};
