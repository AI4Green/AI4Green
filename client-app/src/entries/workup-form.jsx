import "config/i18n";

import { ChakraProvider } from "@chakra-ui/react";
import { WorkupCard } from "components/workups";
import { BackendApiProvider } from "contexts";
import React from "react";
import { createRoot } from "react-dom/client";

/*
 * Mounts a workup card on every `.workup-form-root` on the reaction page (one
 * for the workup, one for the purification). The Jinja template passes the
 * reaction's names and the card's category as data attributes; whether the
 * user can edit comes from the API, which checks it on every write anyway.
 *
 * Known dev-only quirk: StrictMode sets up Chakra v2's modal focus lock twice,
 * so with the Vite dev server a modal briefly takes focus away from its first
 * field (showing its validation message early), focus jumps to the page's
 * "Skip to main content" link when the modal closes, and the modal may not
 * reopen. Production builds don't do the double setup, so none of this happens
 * there.
 */
const CATEGORIES = new Set(["WORKUP", "PURIFICATION"]);

document.querySelectorAll(".workup-form-root").forEach((rootElement) => {
  const { category, reactionId, workgroupName, workbookName } =
    rootElement.dataset;
  if (!CATEGORIES.has(category) || !reactionId) return;

  createRoot(rootElement).render(
    <React.StrictMode>
      <ChakraProvider resetCSS={false}>
        <BackendApiProvider>
          <WorkupCard
            workgroup={workgroupName}
            workbook={workbookName}
            reaction={reactionId}
            category={category}
          />
        </BackendApiProvider>
      </ChakraProvider>
    </React.StrictMode>,
  );
});
