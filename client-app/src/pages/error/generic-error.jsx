import { Container, Text } from "@chakra-ui/react";
// import { TitledAlert } from "components/core/titled-alert";
// import { useTranslation } from "react-i18next";

export const GenericError = () => {
  // const { t } = useTranslation();
  return (
    <Container my={16}>
      {/*<TitledAlert title={t("feedback.error_title")}>*/}
      {/*  <Text>{t("feedback.error")}</Text>*/}
      {/*</TitledAlert>*/}
    </Container>
  );
};
