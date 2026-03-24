import { Container, Text } from "@chakra-ui/react";
// import { TitledAlert } from "components/core/titled-alert";
import { useTranslation } from "react-i18next";

export const Forbidden = () => {
  const { t } = useTranslation();
  return (
    <Container my={16}>
      {/*<TitledAlert title="Forbidden" status="error">*/}
      <Text>{t("feedback.error_403")}</Text>
      {/*</TitledAlert>*/}
    </Container>
  );
};
