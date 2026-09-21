// Generic reusable Busy/Loading fallbacks for Suspense

import { Container } from "@chakra-ui/react";
import { useTranslation } from "react-i18next";

import { LoadingIndicator } from "./loading-indicator";

/**
 * Includes a `Container` so more suitable for Pages
 * @param {object} props
 * @param {string?} [props.tKey] Translation key for loading text
 * @param {object?} [props.containerProps] Props to pass to the container
 * @returns
 */
export const BusyPage = ({
  tKey = "feedback.loading",
  containerProps = {},
}) => {
  const { t } = useTranslation();
  return (
    <Container {...containerProps}>
      <LoadingIndicator verb={t(tKey)} />
    </Container>
  );
};
