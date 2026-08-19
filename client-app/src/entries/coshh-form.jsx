import React, { useState } from "react";
import { createRoot } from "react-dom/client";
import { ChakraProvider } from "@chakra-ui/react";

import { EmbeddedCoshh } from "components/coshh-forms/form";

import { BackendApiProvider } from "contexts";
// import theme from "theme";

const rootElement = document.getElementById("coshh-form-root");

if (rootElement) {
  const reactionId = rootElement.dataset.reactionId;
  const formId = rootElement.dataset.formId || null;

  createRoot(rootElement).render(
    <React.StrictMode>
      <ChakraProvider>
        <BackendApiProvider>
          <EmbeddedCoshh reactionId={reactionId} initialFormId={formId} />
        </BackendApiProvider>
      </ChakraProvider>
    </React.StrictMode>,
  );
}
