import { isFieldTriggered } from "components/section-field";
import { getLatestResponse, initialValues } from "components/section-form";
import { INPUT_TYPES } from "constants";
import DOMPurify from "dompurify";

/*
 * Printable view of a COSHH form instance.
 *
 * Pure rendering only: these components are rendered into the print iframe by
 * export.jsx and styled by coshh-print.css. They use plain HTML elements rather
 * than Chakra, because Chakra's styles are injected into the host page and never
 * reach the iframe.
 */

/**
 * Tags kept by {@link sanitizeHtml}. Formatted text fields hold Quill HTML, so
 * only the inline formats Quill can produce are allowed.
 * @type {string[]}
 */
const ALLOWED_TAGS = [
  "p",
  "br",
  "strong",
  "b",
  "em",
  "i",
  "u",
  "s",
  "sub",
  "sup",
  "span",
];

/**
 * Strips user-supplied HTML down to {@link ALLOWED_TAGS} with DOMPurify. All
 * attributes are removed; disallowed elements are unwrapped so their text is
 * kept, except for elements such as <script> and <style>, which are dropped
 * along with their contents.
 * @param {string} html - Raw HTML from a formatted text field.
 * @returns {string} Sanitised HTML, safe to pass to `dangerouslySetInnerHTML`.
 */
const sanitizeHtml = (html) =>
  DOMPurify.sanitize(html, { ALLOWED_TAGS, ALLOWED_ATTR: [] });

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

    // initialValues always parses file responses into an array
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
 * order. Responses are read with the same helpers as SectionForm, so the
 * printout matches what the user sees on screen.
 * @param {object} props
 * @param {object} props.section - The template section, with its fields and responses.
 * @param {number} props.instanceId - Id of the form instance being printed.
 */
const PrintSection = ({ section, instanceId }) => {
  const fields = section.fields.map((x) => ({
    ...x,
    response: getLatestResponse(x),
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
export const CoshhPrintView = ({ instance, exportedAt }) => (
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
