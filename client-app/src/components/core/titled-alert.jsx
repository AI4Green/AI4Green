import {
  Alert,
  AlertDescription,
  AlertIcon,
  AlertTitle,
  HStack,
  VStack,
} from "@chakra-ui/react";
import { capitalise } from "helpers";

export const TitledAlert = ({ title, status = "info", children, ...p }) => (
  <Alert status={status} {...p}>
    <VStack w="100%" align="start">
      <HStack spacing={0}>
        <AlertIcon boxSize={8} />
        <AlertTitle>{title ?? capitalise(status)}</AlertTitle>
      </HStack>
      {children && <AlertDescription>{children}</AlertDescription>}
    </VStack>
  </Alert>
);
