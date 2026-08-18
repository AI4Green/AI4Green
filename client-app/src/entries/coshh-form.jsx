import React from "react";
import { createRoot } from "react-dom/client";
import { ChakraProvider } from "@chakra-ui/react";

import { CoshhForm } from "components/coshh-forms/form";

import { BackendApiProvider } from "contexts";
// import theme from "theme";

const rootElement = document.getElementById("coshh-form-root");

if (rootElement) {
  const formId = rootElement.dataset.formId;

  createRoot(rootElement).render(
    <React.StrictMode>
      <ChakraProvider>
        <BackendApiProvider>
          <CoshhForm formId={formId} />
        </BackendApiProvider>
      </ChakraProvider>
    </React.StrictMode>,
  );
}
