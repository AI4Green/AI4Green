import { Container, Text } from "@chakra-ui/react";
import { TitledAlert } from "components/core/titled-alert.tsx";
import { useTranslation } from "react-i18next";

export const NotFound = () => {
  const { t } = useTranslation();
  return (
    <Container my={16}>
      <TitledAlert title="404" status="error">
        <Text>{t("feedback.error_404")}</Text>
      </TitledAlert>
    </Container>
  );
};
