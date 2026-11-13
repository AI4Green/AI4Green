import { Box, Flex, Grid, Heading, Icon, Text, VStack } from "@chakra-ui/react";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { useRef, useState } from "react";
import { FaChevronRight } from "react-icons/fa";
import { object, string } from "yup";
import { SelectField, TextField } from "./form-fields.jsx";
import {
  Card,
  EmptyState,
  ErrorMessage,
  INV,
  InvButton,
  MOBILE,
  Muted,
} from "./ui.jsx";
import { describeError } from "./use-inventory-action.js";
import { CAS_PATTERN, formatDate, formatQuantity } from "./utils.jsx";

const FORM_ID = "inv-cheminventory-search";
// the backend's limits, checked here too so a bad search never leaves the page
const MIN_NAME_LENGTH = 2;
const MAX_QUERY_LENGTH = 100;

const SEARCH_TYPES = {
  cas: { label: "CAS number", placeholder: "e.g. 67-64-1" },
  name: { label: "Name", placeholder: "e.g. acetone" },
};

const schema = object({
  type: string().oneOf(Object.keys(SEARCH_TYPES)).required(),
  query: string()
    .trim()
    .required("Enter something to search for")
    .when("type", {
      is: "cas",
      then: (s) => s.matches(CAS_PATTERN, "Enter a CAS number like 67-64-1"),
      otherwise: (s) =>
        s
          .min(MIN_NAME_LENGTH, `Enter at least ${MIN_NAME_LENGTH} characters`)
          .max(
            MAX_QUERY_LENGTH,
            `Enter at most ${MAX_QUERY_LENGTH} characters`,
          ),
    }),
});

/** Formats a count with its noun, e.g. "1 container" or "3 containers". */
const plural = (count, noun) => `${count} ${noun}${count === 1 ? "" : "s"}`;

/**
 * One container in a location tile.
 * @param {object} props
 * @param {object} props.container A container from the search results.
 */
const ContainerRow = ({ container }) => {
  const quantity =
    container.quantity == null
      ? null
      : `${formatQuantity(container.quantity)} ${container.unit}`;
  const details = [
    container.barcode && `Barcode ${container.barcode}`,
    container.supplier,
    container.dateReceived && `Received ${formatDate(container.dateReceived)}`,
  ].filter(Boolean);

  return (
    <Box
      as="li"
      py="10px"
      borderTop="1px solid"
      borderColor={INV.rowLine}
      _first={{ borderTop: "none", pt: 0 }}
    >
      <Flex justify="space-between" gap="12px" align="baseline">
        <Text as="b" m={0}>
          {container.name}
        </Text>
        {quantity && (
          <Text as="span" whiteSpace="nowrap">
            {quantity}
          </Text>
        )}
      </Flex>
      {container.cas && <Muted>CAS {container.cas}</Muted>}
      {details.length > 0 && <Muted>{details.join(" · ")}</Muted>}
    </Box>
  );
};

/**
 * Tile for one location, listing the matching containers held there.
 * @param {object} props
 * @param {object} props.location A location from the search results.
 */
const LocationTile = ({ location }) => {
  const path = location.path?.length ? location.path : [location.name];
  const parents = path.slice(0, -1);

  return (
    <Card
      as="article"
      p="16px"
      borderTop="4px solid"
      borderTopColor={INV.accent}
    >
      {parents.length > 0 && (
        <Flex
          as="ol"
          aria-label="Inside"
          listStyleType="none"
          m={0}
          mb="2px"
          p={0}
          wrap="wrap"
          align="center"
          columnGap="6px"
          color={INV.muted}
          fontSize="12px"
        >
          {parents.map((part, index) => (
            <Flex as="li" key={index} align="center" gap="6px">
              {index > 0 && (
                <Icon as={FaChevronRight} boxSize="8px" aria-hidden="true" />
              )}
              {part}
            </Flex>
          ))}
        </Flex>
      )}
      <Heading
        as="h3"
        m={0}
        fontSize="16px"
        fontWeight={600}
        lineHeight="inherit"
        color={INV.ink}
        title={location.name}
      >
        {path[path.length - 1]}
      </Heading>
      <Muted mb="12px">{plural(location.containers.length, "container")}</Muted>
      <Box as="ul" listStyleType="none" m={0} p={0}>
        {location.containers.map((container) => (
          <ContainerRow key={container.id} container={container} />
        ))}
      </Box>
    </Card>
  );
};

/**
 * The tiles for a search's results, or a message when nothing matched.
 * @param {object} props
 * @param {object} props.results The response from searchChemInventory.
 */
const SearchResults = ({ results }) => {
  const label = SEARCH_TYPES[results.type].label.toLowerCase();

  if (results.locations.length === 0) {
    return (
      <Card>
        <EmptyState>
          No containers in ChemInventory match the {label} “{results.query}”.
        </EmptyState>
      </Card>
    );
  }

  return (
    <VStack
      as="section"
      aria-label="Search results"
      spacing="12px"
      align="stretch"
    >
      <Muted fontSize="14px">
        {plural(results.containerCount, "container")} in{" "}
        {plural(results.locations.length, "location")} match the {label} “
        {results.query}”.
        {results.truncated &&
          " Only the first results are shown, so try a more specific search."}
      </Muted>
      <Grid
        templateColumns="repeat(auto-fill, minmax(280px, 1fr))"
        gap="16px"
        sx={{ [MOBILE]: { gridTemplateColumns: "1fr" } }}
      >
        {results.locations.map((location) => (
          <LocationTile
            key={location.id ?? location.name}
            location={location}
          />
        ))}
      </Grid>
    </VStack>
  );
};

/**
 * Searches a workgroup's ChemInventory by CAS number or name, and shows where
 * the matching containers are, one tile per location. Read-only.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup.
 * @param {object | undefined} props.connection The connection status, from
 *   useChemInventoryConnection; undefined while it loads.
 * @param {() => void} props.onOpenSettings Opens the connection settings.
 */
export default function ChemInventoryView({
  workgroupName,
  connection,
  onOpenSettings,
}) {
  const { inventory: inventoryApi } = useBackendApi();
  const [results, setResults] = useState(null);
  const [error, setError] = useState(null);
  // numbers each search, so a slow earlier search can't replace a later one
  const latestSearch = useRef(0);

  if (!connection) {
    return (
      <Card>
        <EmptyState>Loading…</EmptyState>
      </Card>
    );
  }

  if (!connection.configured) {
    return (
      <Card>
        <EmptyState>
          This workgroup hasn't connected a ChemInventory account.{" "}
          {connection.canManage
            ? "Connect one to search its stock here."
            : "Ask a principal investigator or senior researcher to connect one."}
        </EmptyState>
        {connection.canManage && (
          <Flex justify="center" pb="24px">
            <InvButton primary onClick={onOpenSettings}>
              Connect ChemInventory
            </InvButton>
          </Flex>
        )}
      </Card>
    );
  }

  const handleSearch = async ({ type, query }) => {
    const searchNumber = ++latestSearch.current;
    setError(null);
    try {
      const found = await inventoryApi.searchChemInventory(
        workgroupName,
        type,
        query.trim(),
      );
      if (searchNumber === latestSearch.current) setResults(found);
    } catch (e) {
      const message = await describeError(e);
      if (searchNumber === latestSearch.current) {
        setResults(null);
        setError(message);
      }
    }
  };

  return (
    <>
      <Card p="16px">
        <Formik
          initialValues={{ type: "cas", query: "" }}
          validationSchema={schema}
          onSubmit={handleSearch}
        >
          {({ values, isSubmitting }) => (
            <Flex
              as={Form}
              id={FORM_ID}
              noValidate
              role="search"
              aria-label="Search ChemInventory"
              gap="12px"
              align="flex-start"
              wrap="wrap"
            >
              <Box w="180px" sx={{ [MOBILE]: { w: "100%" } }}>
                <SelectField name="type" label="Search by">
                  {Object.entries(SEARCH_TYPES).map(([value, { label }]) => (
                    <option key={value} value={value}>
                      {label}
                    </option>
                  ))}
                </SelectField>
              </Box>
              <Box flex={1} minW="240px">
                <TextField
                  name="query"
                  label={SEARCH_TYPES[values.type].label}
                  placeholder={SEARCH_TYPES[values.type].placeholder}
                  type="search"
                  autoComplete="off"
                  spellCheck={false}
                  maxLength={MAX_QUERY_LENGTH}
                />
              </Box>
              {/* lines the button up with the inputs, below their labels */}
              <Box pt="23px" sx={{ [MOBILE]: { pt: 0, w: "100%" } }}>
                <InvButton
                  primary
                  type="submit"
                  isLoading={isSubmitting}
                  loadingText="Searching"
                  sx={{ [MOBILE]: { w: "100%" } }}
                >
                  Search
                </InvButton>
              </Box>
            </Flex>
          )}
        </Formik>
        <Muted mt="12px">
          Read-only from ChemInventory
          {connection.inventoryName && (
            <>
              {": "}
              <Text as="b">{connection.inventoryName}</Text>
            </>
          )}
          . Name searches include synonyms.
        </Muted>
      </Card>

      {error && <ErrorMessage mb={0}>{error}</ErrorMessage>}
      {results && <SearchResults results={results} />}
    </>
  );
}
