import { Button, HStack, useToast } from "@chakra-ui/react";
import { fetchKeys as projectFetchKeys } from "api/project";
import { useBackendApi } from "contexts";
import { useState } from "react";
import { flushSync } from "react-dom";
import { createRoot } from "react-dom/client";
import { FaPrint } from "react-icons/fa";
// `?inline` returns the stylesheet as a string without applying it to the app
import printCss from "./coshh-print.css?inline";
import { CoshhPrintView } from "./print-view";

/*
 * Print / PDF export for COSHH form instances.
 *
 * The printable view (print-view.jsx) is built as a standalone HTML document and
 * printed from a hidden iframe, so it never depends on the host page's
 * stylesheets. This keeps it independent of the `@media print` rules in
 * printstyle.css (used by the reaction summary on the reaction constructor page)
 * and means the SPA route and the embedded coshh-form entry behave identically.
 * PDFs come from the browser's "Save as PDF" print destination; the document
 * title becomes the default filename.
 */

/**
 * Builds the print document's title, which browsers use as the default "Save as
 * PDF" filename. Characters outside [\w.-] are replaced with underscores.
 * @param {object} instance - The form instance.
 * @returns {string} A filename-safe title: `COSHH-<reaction code>-<uuid>`.
 */
const exportTitle = (instance) =>
  `COSHH-${instance.reaction_code}-${instance.uuid}`.replace(/[^\w.-]+/g, "_");

/**
 * Builds the empty HTML document loaded into the print frame: the title, the
 * print stylesheet (coshh-print.css, which includes the A4 page setup) and a
 * `#root` container that the print view is rendered into. The title is safe to
 * interpolate because {@link exportTitle} restricts it to [\w.-].
 * @param {object} instance - The form instance.
 * @returns {string} A complete HTML document string.
 */
const renderPrintShell = (instance) =>
  `<!doctype html><html lang="en"><head><meta charset="utf-8">` +
  `<title>${exportTitle(instance)}</title>` +
  `<style>${printCss}</style>` +
  `</head><body><div id="root"></div></body></html>`;

/**
 * Cleanup for the print frame that is currently attached, if any. Kept so the
 * next print can tear the old frame down, because `afterprint` isn't guaranteed
 * to fire in every browser.
 * @type {(() => void) | null}
 */
let cleanupPrintFrame = null;

/**
 * Prints a COSHH form instance from a hidden iframe. The print view is rendered
 * into the frame with React, then the browser's print dialog is opened, where
 * the user can print or choose "Save as PDF". The frame and its React root are
 * removed after printing, or when the next print starts.
 * @param {object} instance - The form instance, as returned by the project API.
 * @returns {Promise<void>} Resolves once the print view has rendered, just
 *   before the print dialog opens; rejects if the print view fails to render.
 */
export const printCoshhInstance = (instance) =>
  new Promise((resolve, reject) => {
    cleanupPrintFrame?.();

    const iframe = document.createElement("iframe");
    iframe.title = "COSHH print preview";
    iframe.setAttribute("aria-hidden", "true");
    iframe.tabIndex = -1;
    // No `allow-scripts`: nothing inside the frame can run, even if user HTML
    // got past the sanitiser. `allow-same-origin` lets this page render into
    // the frame and `allow-modals` lets it open the print dialog.
    iframe.setAttribute("sandbox", "allow-same-origin allow-modals");
    iframe.style.cssText =
      "position:fixed;right:0;bottom:0;width:0;height:0;border:0;visibility:hidden;";

    // Declared here, assigned in onload, so cleanup can unmount it if it exists.
    let root = null;
    const cleanup = () => {
      root?.unmount();
      iframe.remove();
      // Only clear the shared slot if it still points at this frame, so a late
      // cleanup from an old frame can't forget a newer one.
      if (cleanupPrintFrame === cleanup) cleanupPrintFrame = null;
    };
    cleanupPrintFrame = cleanup;

    iframe.onload = () => {
      try {
        const frameWindow = iframe.contentWindow;
        root = createRoot(frameWindow.document.getElementById("root"));
        // render synchronously so the DOM is complete before print() snapshots it
        flushSync(() => {
          root.render(
            <CoshhPrintView instance={instance} exportedAt={new Date()} />,
          );
        });
        // Defer cleanup so the frame isn't removed while its own print() call is
        // still on the stack, which can stop the rest of this handler running.
        frameWindow.addEventListener("afterprint", () => setTimeout(cleanup), {
          once: true,
        });
        // Resolve before printing: print() blocks until the dialog closes in
        // most browsers, and callers shouldn't depend on it ever returning.
        resolve();
        frameWindow.focus();
        frameWindow.print();
      } catch (e) {
        cleanup();
        reject(e);
      }
    };

    iframe.srcdoc = renderPrintShell(instance);
    document.body.appendChild(iframe);
  });

/**
 * "Print to PDF" button for a COSHH form. On click it fetches the latest copy of
 * the form (the SWR cache doesn't see autosaved changes) and passes it to
 * {@link printCoshhInstance}. Shows a loading state while working and an error
 * toast if the export fails.
 * @param {object} props
 * @param {number} props.formId - Id of the form instance to export.
 */
export const CoshhExportActions = ({ formId }) => {
  const { apiFetcher } = useBackendApi();
  const toast = useToast();
  const [busy, setBusy] = useState(false);

  const handlePrint = async () => {
    try {
      setBusy(true);
      // fetch fresh, the SWR cache doesn't see autosaved changes
      const instance = await apiFetcher(projectFetchKeys.get(formId));
      await printCoshhInstance(instance);
    } catch (e) {
      console.error(e);
      toast({
        title: "Failed to export COSHH form",
        status: "error",
        position: "top",
        isClosable: true,
      });
    } finally {
      setBusy(false);
    }
  };

  return (
    <HStack justify="flex-end" spacing={2} mt={4} mb={4} pr={4}>
      <Button
        size="sm"
        colorScheme="green"
        leftIcon={<FaPrint />}
        onClick={handlePrint}
        isLoading={busy}
      >
        Print to PDF
      </Button>
    </HStack>
  );
};
