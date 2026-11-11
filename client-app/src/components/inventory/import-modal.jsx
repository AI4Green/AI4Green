import {
  Box,
  FormControl,
  FormHelperText,
  FormLabel,
  Input,
  Select,
  Table,
  Tbody,
  Text,
  Thead,
  Tr,
  VStack,
  useToast,
} from "@chakra-ui/react";
import { useWorkgroupInventory } from "api";
import { GLOBAL_PARAMETERS } from "constants";
import { useBackendApi } from "contexts";
import { useState } from "react";
import {
  Badge,
  ErrorMessage,
  INV,
  InvButton,
  InvModal,
  InvTd,
  InvTh,
  Muted,
  controlProps,
  tableProps,
} from "./ui.jsx";
import { describeError } from "./use-inventory-action.js";

// matches MAX_FILE_BYTES in services/inventory_import.py; the server checks too
const MAX_FILE_BYTES = 1024 * 1024;

/** The CSV columns, in template order, with what users need to know. */
const COLUMNS = [
  {
    name: "CAS number",
    required: false,
    help: "Matched to the compound database. Leave blank if there isn't one.",
  },
  {
    name: "Name",
    required: true,
    help: "Used to match a compound with no CAS.",
  },
  {
    name: "SMILES",
    required: false,
    help: "Required for compounds that aren't in the database. Also used to match compounds by structure.",
  },
  {
    name: "Hazard codes",
    required: false,
    help: "For new compounds only, e.g. H225-H319. Statements are added automatically.",
  },
  {
    name: "Location",
    required: true,
    help: "Matched by name; created if it doesn't exist.",
  },
  {
    name: "Location description",
    required: false,
    help: "Used when the import creates the location.",
  },
  { name: "Quantity", required: true, help: "A number, e.g. 500." },
  { name: "Unit", required: true, help: "mg, g, kg, µL, mL or L." },
  { name: "Supplier", required: false },
  { name: "Lot number", required: false },
  { name: "Purity", required: false },
  {
    name: "Barcode",
    required: false,
    help: "Rows whose barcode is already in the inventory are skipped, so the same file can be uploaded again safely.",
  },
  { name: "Date received", required: false, help: "YYYY-MM-DD" },
  { name: "Expiry date", required: false, help: "YYYY-MM-DD" },
];

/** Downloads a CSV with just the header row. */
const downloadTemplate = () => {
  const csv = `${COLUMNS.map((column) => column.name).join(",")}\r\n`;
  const url = URL.createObjectURL(new Blob([csv], { type: "text/csv" }));
  const link = document.createElement("a");
  link.href = url;
  link.download = "ai4green-inventory-import-template.csv";
  link.click();
  // revoke after the click has started the download
  setTimeout(() => URL.revokeObjectURL(url), 0);
};

/**
 * Checks the chosen file before it's sent, so obvious mistakes are caught
 * straight away. The server validates it properly either way.
 * @param {File} file The chosen file.
 * @returns {string | null} What's wrong, or null.
 */
const checkFile = (file) => {
  if (!file.name.toLowerCase().endsWith(".csv")) {
    return "Choose a .csv file. In Excel, use Save As and pick CSV UTF-8.";
  }
  if (file.size > MAX_FILE_BYTES) {
    return "The file is over 1 MB. Split it into smaller files.";
  }
  return null;
};

/**
 * Reads the per-row results from a failed import, which the backend sends
 * when any row is invalid.
 * @param {Error} error The error thrown by ky.
 * @returns {Promise<object | null>} The results, or null if there aren't any.
 */
const readImportResult = async (error) => {
  if (error?.response?.status !== 400) return null;
  try {
    const body = await error.response.clone().json();
    return Array.isArray(body?.rows) ? body : null;
  } catch {
    return null;
  }
};

/**
 * The list of columns users need to fill in.
 */
const ColumnHelp = () => (
  <Box as="ul" listStyleType="none" m={0} fontSize="13px">
    {COLUMNS.map((column) => (
      <Box as="li" key={column.name} mb="4px">
        <Text as="b">{column.name}</Text>
        {column.required ? " (required)" : ""}
        {column.help && (
          <Text as="span" color={INV.muted}>
            {" "}
            {column.help}
          </Text>
        )}
      </Box>
    ))}
  </Box>
);

/**
 * What the import would do, with every row that has an error or warning.
 * @param {object} props
 * @param {object} props.result The preview from the backend.
 */
const PreviewResult = ({ result }) => {
  const { summary, rows } = result;
  // skipped duplicates are only counted, so re-uploading a large file
  // doesn't list every row
  const flagged = rows.filter(
    (row) => row.error || (!row.duplicate && row.warnings.length > 0),
  );

  return (
    <VStack align="stretch" spacing="12px">
      <Box role="status">
        {summary.errorCount > 0 ? (
          <ErrorMessage mb={0}>
            {summary.errorCount} of {summary.rowCount} rows have problems. Fix
            them in the file and check it again. Nothing has been imported.
          </ErrorMessage>
        ) : summary.newItemCount === 0 ? (
          <Text>
            Everything in this file is already in the inventory, so there's
            nothing to import.
          </Text>
        ) : (
          <Text>
            Ready to import <Text as="b">{summary.newItemCount}</Text>{" "}
            containers, adding {summary.newLocationCount} new locations and{" "}
            {summary.newCompoundCount} new compounds.
          </Text>
        )}
        {summary.errorCount === 0 && summary.duplicateCount > 0 && (
          <Muted>
            {summary.duplicateCount} rows are already in the inventory and will
            be skipped.
          </Muted>
        )}
      </Box>
      {flagged.length > 0 && (
        <Table {...tableProps}>
          <Thead>
            <Tr>
              <InvTh>Row</InvTh>
              <InvTh>Compound</InvTh>
              <InvTh>Problem</InvTh>
            </Tr>
          </Thead>
          <Tbody>
            {flagged.map((row) => (
              <Tr key={row.row}>
                <InvTd>{row.row}</InvTd>
                <InvTd>{row.compound ?? row.name}</InvTd>
                <InvTd>
                  {row.error && <Badge tone="exp">{row.error}</Badge>}
                  {row.warnings.map((warning) => (
                    <Muted key={warning}>{warning}</Muted>
                  ))}
                </InvTd>
              </Tr>
            ))}
          </Tbody>
        </Table>
      )}
    </VStack>
  );
};

/**
 * Imports containers from a CSV file. The file is checked first and the user
 * sees what will happen; nothing is saved until they confirm, and then only
 * if every row is valid.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup to import into.
 * @param {object[]} props.workbooks The user's workbooks, as `{ id, name }`,
 *   for compounds that aren't in the database.
 * @param {() => void} props.onClose Closes the modal.
 */
export const ImportModal = ({ workgroupName, workbooks, onClose }) => {
  const { inventory } = useBackendApi();
  const { mutate } = useWorkgroupInventory(workgroupName);
  const toast = useToast();
  const [file, setFile] = useState(null);
  const [workbookId, setWorkbookId] = useState(
    workbooks.length === 1 ? workbooks[0].id : null,
  );
  const [result, setResult] = useState(null);
  const [error, setError] = useState(null);
  const [isLoading, setIsLoading] = useState(false);

  const canImport =
    result &&
    result.summary.errorCount === 0 &&
    result.summary.newItemCount > 0;

  const resetResult = () => {
    setResult(null);
    setError(null);
  };

  const handleFileChange = (event) => {
    resetResult();
    const chosen = event.target.files?.[0] ?? null;
    const problem = chosen && checkFile(chosen);
    setFile(problem ? null : chosen);
    if (problem) setError(problem);
  };

  const handleWorkbookChange = (event) => {
    resetResult();
    setWorkbookId(event.target.value ? Number(event.target.value) : null);
  };

  const handleCheck = async () => {
    resetResult();
    if (!file) {
      setError("Choose a CSV file to import.");
      return;
    }
    setIsLoading(true);
    try {
      setResult(await inventory.previewImport(workgroupName, file, workbookId));
    } catch (e) {
      console.error(e);
      setError(await describeError(e));
    } finally {
      setIsLoading(false);
    }
  };

  const handleImport = async () => {
    setIsLoading(true);
    setError(null);
    try {
      const imported = await inventory.importInventory(
        workgroupName,
        file,
        workbookId,
      );
      await mutate();
      toast({
        title: `Imported ${imported.summary.newItemCount} containers`,
        description:
          imported.summary.duplicateCount > 0
            ? `Skipped ${imported.summary.duplicateCount} already in the inventory.`
            : undefined,
        status: "success",
        duration: GLOBAL_PARAMETERS.ToastDuration,
        isClosable: true,
        position: "top",
      });
      onClose();
    } catch (e) {
      console.error(e);
      // the stock may have changed since the preview, so show the new results
      const rowResult = await readImportResult(e);
      if (rowResult) setResult(rowResult);
      else setError(await describeError(e));
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <InvModal
      isOpen
      onClose={onClose}
      title="Import stock from CSV"
      size="2xl"
      actionLabel={canImport ? "Import" : "Check file"}
      onAction={canImport ? handleImport : handleCheck}
      isLoading={isLoading}
    >
      {error && <ErrorMessage>{error}</ErrorMessage>}
      <VStack align="stretch" spacing="16px">
        <Text>
          Each row is one container. Compounds are matched by CAS number, then
          by structure (SMILES), then by name if there's no CAS number.
        </Text>
        <ColumnHelp />
        <Box>
          <InvButton onClick={downloadTemplate}>Download template</InvButton>
        </Box>
        <FormControl>
          <FormLabel fontSize="14px">Workbook for new compounds</FormLabel>
          <Select
            {...controlProps}
            placeholder="None"
            value={workbookId ?? ""}
            onChange={handleWorkbookChange}
          >
            {workbooks.map((workbook) => (
              <option key={workbook.id} value={workbook.id}>
                {workbook.name}
              </option>
            ))}
          </Select>
          <FormHelperText fontSize="12px">
            Compounds that aren't in the database are added here as novel
            compounds, using their SMILES.
          </FormHelperText>
        </FormControl>
        <FormControl isRequired>
          <FormLabel fontSize="14px">CSV file</FormLabel>
          <Input
            type="file"
            accept=".csv,text/csv"
            onChange={handleFileChange}
            p="8px"
            h="auto"
            borderColor={INV.lineStrong}
          />
        </FormControl>
        {result && <PreviewResult result={result} />}
      </VStack>
    </InvModal>
  );
};
