// import "config/i18n";

import { ChakraProvider, Flex, Spinner } from "@chakra-ui/react";
import { StrictMode, Suspense } from "react";
import { createRoot } from "react-dom/client";
import { BrowserRouter } from "react-router-dom";

// import { Root } from "./routes/root";
import { theme } from "./themes";
import App from "./App.tsx";

createRoot(document.getElementById("root")).render(
  <StrictMode>
    <ChakraProvider value={theme}>
      <BrowserRouter>
        <Suspense
          fallback={
            <Flex justify="center" w="100%" my={16}>
              <Spinner boxSize={16} />
            </Flex>
          }
        >
          <App />
        </Suspense>
      </BrowserRouter>
    </ChakraProvider>
  </StrictMode>,
);
