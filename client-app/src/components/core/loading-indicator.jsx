import { Flex, Spinner, Text } from "@chakra-ui/react";
import { useTranslation } from "react-i18next";

export const LoadingIndicator = ({ verb, noun, layoutProps, textProps }) => {
  const { t } = useTranslation();

  verb = verb ?? t("feedback.loading");

  return (
    <Flex justify="center" {...layoutProps}>
      <Flex justify="space-evenly" p={5} {...textProps}>
        <Flex mr={2}>
          <Spinner />
        </Flex>
        <Text as="i">
          {verb} {noun || null} ...
        </Text>
      </Flex>
    </Flex>
  );
};
