import {
  Alert,
  AlertIcon,
  Button,
  Grid,
  GridItem,
  HStack,
  Icon,
  Stack,
  Text,
  VStack,
} from "@chakra-ui/react";
import { Badge } from "components/core/Badge.js";
import { Footer } from "components/core/footer.js";
import { Sidebar } from "components/core/nav";
import { useLocationStateToast } from "helpers/hooks";
import { FaExclamationTriangle } from "react-icons/fa";
import { IoIosAddCircleOutline } from "react-icons/io";
import { Link, Outlet } from "react-router-dom";

export const DefaultLayout = ({
  toastDefaults = { position: "top" },
  children,
}) => {
  useLocationStateToast(toastDefaults);

  return (
    <Grid templateRows="1fr auto" minHeight="100vh" fontWeight="light">
      <Sidebar>{children ? children : <Outlet />}</Sidebar>
      <GridItem>{<Footer />}</GridItem>
    </Grid>
  );
};

export const DefaultContentLayout = ({ children }) => (
  <Stack align="center" minW="full" pb={8}>
    <VStack
      p={4}
      w={{ base: "full", xl: "90%", "2xl": "70%" }}
      spacing={4}
      align="stretch"
      borderWidth={0.5}
      borderRadius={4}
    >
      {children}
    </VStack>
  </Stack>
);

export const DefaultContentHeader = ({ icon, header }) => (
  <Badge
    label={header}
    leftIcon={icon}
    colorScheme="brand"
    variant="outline"
    fontSize="xs"
  />
);

export const NewButton = ({
  icon = IoIosAddCircleOutline,
  label = "New",
  onClick,
  to,
}) => (
  <Button
    as={onClick ? undefined : Link}
    onClick={onClick}
    to={to}
    colorScheme="gray"
    leftIcon={<Icon as={icon} fontSize={14} color="green" />}
    size="xs"
    borderRadius={10}
    variant="outline"
    px={4}
    py={3.5}
  >
    <Text
      fontSize={{ base: "xs", md: "sm" }}
      fontWeight="normal"
      color="gray.700"
    >
      {label}
    </Text>
  </Button>
);

export const ConfirmationModal = ({
  iconProps = {
    Icon: FaExclamationTriangle,
    color: "red.500",
    fontSize: "4xl",
  },
  feedback,
  content,
}) => {
  return (
    <HStack spacing={4}>
      <Icon
        as={iconProps.Icon}
        color={iconProps.color}
        fontSize={iconProps.fontSize}
      />

      <VStack align="flex-start" flex={1}>
        {feedback && (
          <Alert status={feedback.status}>
            <AlertIcon />
            {feedback.message}
          </Alert>
        )}
        <Text fontSize="sm" fontWeight="light">
          {content.description}
        </Text>

        <VStack
          align="flex-start"
          borderWidth={1}
          borderRadius={7}
          p={2}
          w="full"
          spacing={2}
        >
          <Text fontWeight="semibold" fontSize="sm">
            {content.value}
          </Text>
          {content.tags?.length > 0 && (
            <HStack>
              {content.tags.map((tag) => (
                <Badge
                  key={tag.label}
                  colorScheme={tag.colorScheme}
                  label={tag.label}
                  leftIcon={tag.leftIcon}
                  variant="outline"
                  fontSize="xxs"
                />
              ))}
            </HStack>
          )}
        </VStack>
      </VStack>
    </HStack>
  );
};
