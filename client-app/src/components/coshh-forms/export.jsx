import { Button, HStack, useToast } from "@chakra-ui/react";
import { fetchKeys as projectFetchKeys } from "api/project";
import { isFieldTriggered } from "components/section-field";
import { initialValues } from "components/section-form";
import { INPUT_TYPES } from "constants";
import { useBackendApi } from "contexts";
import { useState } from "react";
import { flushSync } from "react-dom";
import { createRoot } from "react-dom/client";
import { FaPrint } from "react-icons/fa";

/*
 * Print / PDF export for COSHH form instances.
 *
 * The printable view is built as a standalone HTML document and printed from a
 * hidden iframe, so it never depends on the host page's stylesheets. This keeps it
 * independent of the `@media print` rules in printstyle.css (used by the reaction
 * summary on the reaction constructor page) and means the SPA route and the
 * embedded coshh-form entry behave identically. PDFs come from the browser's
 * "Save as PDF" print destination; the document title becomes the default filename.
 */

/**
 * Stylesheet for the printable view, injected into the print frame's <head>.
 * Every rule is scoped to `.coshh-print` and sized in points for A4 output.
 * @type {string}
 */
const PRINT_CSS = `
.coshh-print { font-family: Arial, Helvetica, sans-serif; font-size: 11pt; color: #000; background: #fff; line-height: 1.4; -webkit-print-color-adjust: exact; print-color-adjust: exact; }
.coshh-print *, .coshh-print *::before, .coshh-print *::after { box-sizing: border-box; }
.coshh-print h1 { font-size: 18pt; margin: 0 0 4pt; }
.coshh-print .coshh-section { border: 1pt solid #888; border-radius: 4pt; padding: 0 8pt 8pt; margin: 12pt 0 0; }
.coshh-print h2 { font-size: 13pt; margin: 0 -8pt 8pt; padding: 5pt 8pt; background: #f0f0f0; border-bottom: 1pt solid #888; border-radius: 4pt 4pt 0 0; page-break-after: avoid; break-after: avoid; }
.coshh-print h3 { font-size: 11pt; margin: 10pt 0 4pt; text-decoration: underline; page-break-after: avoid; break-after: avoid; }
.coshh-print .coshh-meta { width: 100%; border-collapse: collapse; margin-bottom: 8pt; }
.coshh-print .coshh-meta th { text-align: left; width: 28%; font-weight: bold; }
.coshh-print .coshh-meta th, .coshh-print .coshh-meta td { padding: 2pt 4pt; border: 0.75pt solid #999; }
.coshh-print .coshh-field { margin: 0 0 8pt; page-break-inside: avoid; break-inside: avoid; }
.coshh-print .coshh-label { font-weight: bold; margin-bottom: 2pt; }
.coshh-print .coshh-value { border: 0.75pt solid #999; border-radius: 3pt; padding: 4pt 6pt; min-height: 18pt; white-space: pre-wrap; overflow-wrap: anywhere; word-break: break-word; }
.coshh-print .coshh-value p { margin: 0; }
.coshh-print .coshh-value table.coshh-table { margin-top: 0; }
.coshh-print .coshh-content { border: 0.75pt solid #ccc; border-radius: 3pt; padding: 4pt 6pt; background: #fafafa; white-space: pre-wrap; overflow-wrap: anywhere; }
.coshh-print .coshh-empty { color: #777; font-style: italic; }
.coshh-print ul { margin: 2pt 0; padding-left: 16pt; }
.coshh-print table.coshh-table { width: 100%; border-collapse: collapse; margin-top: 2pt; font-size: 9.5pt; }
.coshh-print table.coshh-table th, .coshh-print table.coshh-table td { border: 0.75pt solid #999; padding: 2pt 4pt; text-align: left; vertical-align: top; white-space: pre-wrap; overflow-wrap: anywhere; }
.coshh-print table.coshh-table th { background: #eee; }
.coshh-print tr { page-break-inside: avoid; break-inside: avoid; }
.coshh-print .coshh-footer { margin-top: 16pt; font-size: 8.5pt; color: #555; }
`;

/**
 * Tags kept by {@link sanitizeHtml}. Formatted text fields hold Quill HTML, so
 * only the inline formats Quill can produce are allowed.
 * @type {Set<string>}
 */
const ALLOWED_TAGS = new Set([
  "P",
  "BR",
  "STRONG",
  "B",
  "EM",
  "I",
  "U",
  "S",
  "SUB",
  "SUP",
  "SPAN",
]);

/**
 * Strips user-supplied HTML down to {@link ALLOWED_TAGS}. Disallowed elements
 * are unwrapped (their text is kept), all attributes are removed, and comments
 * and other non-element nodes are dropped.
 * @param {string} html - Raw HTML from a formatted text field.
 * @returns {string} Sanitised HTML, safe to pass to `dangerouslySetInnerHTML`.
 */
const sanitizeHtml = (html) => {
  const doc = new DOMParser().parseFromString(html, "text/html");

  const clean = (node) => {
    [...node.childNodes].forEach((child) => {
      if (child.nodeType === Node.TEXT_NODE) return;
      if (child.nodeType !== Node.ELEMENT_NODE) {
        child.remove();
        return;
      }
      clean(child);
      if (ALLOWED_TAGS.has(child.tagName)) {
        [...child.attributes].forEach((attr) =>
          child.removeAttribute(attr.name),
        );
      } else {
        child.replaceWith(...child.childNodes);
      }
    });
  };

  clean(doc.body);
  return doc.body.innerHTML;
};

/**
 * Whether a response value counts as unanswered.
 * @param {*} value - A field response value.
 * @returns {boolean} True for null, undefined, an empty string or an empty array.
 */
const isEmpty = (value) =>
  value === null ||
  value === undefined ||
  value === "" ||
  (Array.isArray(value) && value.length === 0);

/**
 * Turns an object key into a readable column heading,
 * e.g. `"hazardCode"` or `"hazard_code"` becomes `"Hazard code"`.
 * @param {string} key - A camelCase, snake_case or kebab-case key.
 * @returns {string} The key as a sentence-case label.
 */
const humanise = (key) =>
  String(key)
    .replace(/([a-z])([A-Z])/g, "$1 $2")
    .replace(/[_-]+/g, " ")
    .replace(/^\w/, (c) => c.toUpperCase());

/**
 * Formats a single value for display. Booleans become "Yes"/"No"; objects use
 * their `name`, `label` or `display_name`, falling back to JSON.
 * @param {*} value - The value to format.
 * @returns {string} Display text.
 */
const formatScalar = (value) => {
  if (typeof value === "boolean") return value ? "Yes" : "No";
  if (value && typeof value === "object") {
    return (
      value.name ?? value.label ?? value.display_name ?? JSON.stringify(value)
    );
  }
  return String(value);
};

/**
 * Fallback renderer for table-style and other structured responses.
 * - A plain object renders as a two-column key/value table.
 * - An array of objects renders as a table with one column per key, skipping
 *   internal identifiers and binary/image payloads.
 * - Any other array renders as a bulleted list.
 * - A scalar renders as text.
 * @param {object} props
 * @param {*} props.value - The response value.
 */
const GenericValue = ({ value }) => {
  if (!Array.isArray(value)) {
    if (value && typeof value === "object") {
      return (
        <table className="coshh-table">
          <tbody>
            {Object.entries(value).map(([key, v]) => (
              <tr key={key}>
                <th>{humanise(key)}</th>
                <td>{isEmpty(v) ? "" : formatScalar(v)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      );
    }
    return formatScalar(value);
  }

  const rows = value.filter((row) => !isEmpty(row));
  const isObjectRows = rows.every(
    (row) => row && typeof row === "object" && !Array.isArray(row),
  );

  if (!isObjectRows) {
    return (
      <ul>
        {rows.map((row, i) => (
          <li key={i}>{formatScalar(row)}</li>
        ))}
      </ul>
    );
  }

  const columns = [...new Set(rows.flatMap((row) => Object.keys(row)))].filter(
    (key) => !/^(id|uuid|instanceId|isNew|file|image)$/i.test(key),
  );

  return (
    <table className="coshh-table">
      <thead>
        <tr>
          {columns.map((key) => (
            <th key={key}>{humanise(key)}</th>
          ))}
        </tr>
      </thead>
      <tbody>
        {rows.map((row, i) => (
          <tr key={i}>
            {columns.map((key) => (
              <td key={key}>
                {isEmpty(row[key]) ? "" : formatScalar(row[key])}
              </td>
            ))}
          </tr>
        ))}
      </tbody>
    </table>
  );
};

/**
 * Renders a field's response according to its input type, showing
 * "Not provided" when the field is unanswered. Unknown types fall back to
 * {@link GenericValue}.
 * @param {object} props
 * @param {object} props.field - The template field, including its `inputType`.
 * @param {*} props.value - The field's current response value.
 */
const FieldValue = ({ field, value }) => {
  if (isEmpty(value)) return <span className="coshh-empty">Not provided</span>;

  switch (field.inputType.title.toLowerCase()) {
    case INPUT_TYPES.Radio:
    case INPUT_TYPES.Multiple:
      return Array.isArray(value)
        ? value.map((option) => option.name ?? option).join(", ")
        : formatScalar(value);

    case INPUT_TYPES.File:
    case INPUT_TYPES.ImageFile:
      return (
        <ul>
          {value.map((file) => (
            <li key={file.uuid ?? file.name}>{file.name}</li>
          ))}
        </ul>
      );

    case INPUT_TYPES.DateAndTime: {
      const date = new Date(value);
      return Number.isNaN(date.getTime()) ? value : date.toLocaleString();
    }

    case INPUT_TYPES.FormattedTextInput:
      return (
        <div
          dangerouslySetInnerHTML={{ __html: sanitizeHtml(String(value)) }}
        />
      );

    case INPUT_TYPES.Text:
    case INPUT_TYPES.Description:
    case INPUT_TYPES.Number:
      return String(value);

    default:
      return <GenericValue value={value} />;
  }
};

/**
 * Renders one field of a section: a sub-heading for header fields, the fixed
 * text for content fields, or a labelled response for everything else. If the
 * field's trigger condition is met, its target field is rendered after it,
 * mirroring SectionField.
 * @param {object} props
 * @param {object} props.field - The field to render.
 * @param {object[]} props.fields - All fields in the section, used to find trigger targets.
 * @param {Object<string, *>} props.values - Response values keyed by field id.
 */
const PrintField = ({ field, fields, values }) => {
  const type = field.inputType.title.toLowerCase();

  let body;
  if (type === INPUT_TYPES.Header) {
    body = <h3>{field.name}</h3>;
  } else if (type === INPUT_TYPES.Content) {
    body = (
      <div className="coshh-field">
        <div className="coshh-label">{field.name}</div>
        <div className="coshh-content">{field.defaultResponse}</div>
      </div>
    );
  } else {
    body = (
      <div className="coshh-field">
        <div className="coshh-label">{field.name}</div>
        <div className="coshh-value">
          <FieldValue field={field} value={values[field.id]} />
        </div>
      </div>
    );
  }

  const target =
    field.triggerField &&
    isFieldTriggered(
      field.inputType.name,
      field.triggerField.value,
      values[field.id] ?? [],
    ) &&
    fields.find((x) => x.id === field.triggerField.id);

  return (
    <>
      {body}
      {target && <PrintField field={target} fields={fields} values={values} />}
    </>
  );
};

/**
 * Renders one form section as a bordered block with its visible fields in sort
 * order. Responses are mapped the same way as SectionForm so the printout
 * matches what the user sees on screen.
 * @param {object} props
 * @param {object} props.section - The template section, with its fields and responses.
 * @param {number} props.instanceId - Id of the form instance being printed.
 */
const PrintSection = ({ section, instanceId }) => {
  const fields = section.fields.map((x) => ({
    ...x,
    response: {
      id: x.fieldResponse?.[0]?.id || null,
      value: x.fieldResponse?.[0]?.fieldResponseValues?.at(-1)?.value || null,
    },
  }));
  const values = initialValues(fields, instanceId, section.id);

  return (
    <section className="coshh-section">
      <h2>{section.name}</h2>
      {[...fields]
        .sort((a, b) => a.sortOrder - b.sortOrder)
        .map(
          (field) =>
            !field.hidden && (
              <PrintField
                key={field.id}
                field={field}
                fields={fields}
                values={values}
              />
            ),
        )}
    </section>
  );
};

/**
 * The complete printable COSHH assessment: a details table (template, reaction,
 * workgroup/workbook, owner, form id), every section in sort order, and an
 * export timestamp footer.
 * @param {object} props
 * @param {object} props.instance - The form instance, as returned by the project API.
 * @param {Date} props.exportedAt - When the export was requested, shown in the footer.
 */
const CoshhPrintView = ({ instance, exportedAt }) => (
  <div className="coshh-print">
    <h1>COSHH Assessment</h1>
    <table className="coshh-meta">
      <tbody>
        <tr>
          <th>Template</th>
          <td>{instance.template_name}</td>
        </tr>
        <tr>
          <th>Reaction</th>
          <td>{instance.reaction_code}</td>
        </tr>
        <tr>
          <th>Workgroup / Workbook</th>
          <td>
            {instance.workgroup} / {instance.workbook}
          </td>
        </tr>
        <tr>
          <th>Owner</th>
          <td>{instance.owner_email}</td>
        </tr>
        <tr>
          <th>Form ID</th>
          <td>{instance.uuid}</td>
        </tr>
      </tbody>
    </table>
    {[...(instance.sections ?? [])]
      .sort((a, b) => a.sortOrder - b.sortOrder)
      .map((section) => (
        <PrintSection
          key={section.id}
          section={section}
          instanceId={instance.id}
        />
      ))}
    <div className="coshh-footer">
      Exported from AI4Green on {exportedAt.toLocaleString()}
    </div>
  </div>
);

/**
 * Builds the print document's title, which browsers use as the default "Save as
 * PDF" filename. Characters outside [\w.-] are replaced with underscores.
 * @param {object} instance - The form instance.
 * @returns {string} A filename-safe title: `COSHH-<reaction code>-<uuid>`.
 */
const exportTitle = (instance) =>
  `COSHH-${instance.reaction_code}-${instance.uuid}`.replace(/[^\w.-]+/g, "_");

/**
 * Builds the empty HTML document loaded into the print frame: the title, A4 page
 * setup, {@link PRINT_CSS} and a `#root` container that the print view is
 * rendered into. The title is safe to interpolate because {@link exportTitle}
 * restricts it to [\w.-].
 * @param {object} instance - The form instance.
 * @returns {string} A complete HTML document string.
 */
const renderPrintShell = (instance) =>
  `<!doctype html><html lang="en"><head><meta charset="utf-8">` +
  `<title>${exportTitle(instance)}</title>` +
  `<style>@page { size: A4; margin: 15mm; } body { margin: 0; }${PRINT_CSS}</style>` +
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
 * @returns {Promise<void>} Resolves once the print dialog has been opened;
 *   rejects if the print view fails to render.
 */
export const printCoshhInstance = (instance) =>
  new Promise((resolve, reject) => {
    cleanupPrintFrame?.();

    const iframe = document.createElement("iframe");
    iframe.title = "COSHH print preview";
    iframe.setAttribute("aria-hidden", "true");
    iframe.tabIndex = -1;
    iframe.style.cssText =
      "position:fixed;right:0;bottom:0;width:0;height:0;border:0;visibility:hidden;";

    let root = null;
    const cleanup = () => {
      root?.unmount();
      iframe.remove();
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
        frameWindow.addEventListener("afterprint", cleanup, { once: true });
        frameWindow.focus();
        frameWindow.print();
        resolve();
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
