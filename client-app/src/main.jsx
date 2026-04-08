import "config/i18n";

import { ChakraProvider, Flex, Spinner } from "@chakra-ui/react";
import { StrictMode, Suspense } from "react";
import { createRoot } from "react-dom/client";
import { BrowserRouter } from "react-router-dom";
import {
  BackendApiProvider,
  BackendConfigProvider,
  UserProvider,
} from "./contexts";
import { Root } from "./routes/root";
import { theme } from "./themes";

createRoot(document.getElementById("root")).render(
  <StrictMode>
    <ChakraProvider value={theme}>
      <BrowserRouter basename="/spa">
        <Suspense
          fallback={
            <Flex justify="center" w="100%" my={16}>
              <Spinner boxSize={16} />
            </Flex>
          }
        >
          <BackendApiProvider>
            <UserProvider>
              <BackendConfigProvider>
                <Root />
              </BackendConfigProvider>
            </UserProvider>
          </BackendApiProvider>
        </Suspense>
      </BrowserRouter>
    </ChakraProvider>
  </StrictMode>,
);
