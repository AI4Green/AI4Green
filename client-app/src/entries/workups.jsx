import React, { useState } from "react";
import { createRoot } from "react-dom/client";
import { ChakraProvider } from "@chakra-ui/react";

import { EmbeddedCoshh } from "components/coshh-forms/form";

import { BackendApiProvider } from "contexts";
// import theme from "theme";

const rootElement = document.getElementById("workup-route");

if (rootElement) {
  const reactionId = rootElement.dataset.reactionId;
  const workgroupName = rootElement.dataset.workgroupName;
  const workbookName = rootElement.dataset.workbookName;
  const instanceId = rootElement.dataset.instanceId || null;

  createRoot(rootElement).render(
    <React.StrictMode>
      <ChakraProvider resetCSS={false}>
        <BackendApiProvider>
          <EmbeddedWorkup
            reactionId={reactionId}
            initialFormId={instanceId}
            workgroupName={workgroupName}
            workbookName={workbookName}
          />
        </BackendApiProvider>
      </ChakraProvider>
    </React.StrictMode>,
  );
}
