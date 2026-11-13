import {
  Box,
  Checkbox,
  CheckboxGroup,
  Flex,
  FormControl,
  FormLabel,
  Grid,
  Heading,
  Stack,
  VStack,
} from "@chakra-ui/react";
import { STRUCTURE_SEARCH_ROLES, useSearchWorkgroups } from "api/search";
import { StructureSketcher } from "components/core/structure-sketcher/structure-sketcher.jsx";
import {
  ROLE_LABELS,
  STRUCTURE_MODE_LABELS,
} from "components/core/structure-sketcher/structure-search-labels.js";
import { SelectField, TextField } from "components/inventory/form-fields.jsx";
import { SegmentButton } from "components/inventory/inventory-filters.jsx";
import {
  Card,
  EmptyState,
  ErrorMessage,
  INV,
  InvButton,
  InvPage,
  Muted,
} from "components/inventory/ui.jsx";
import { describeError } from "components/inventory/use-inventory-action.js";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { useRef, useState } from "react";
import { array, object, string } from "yup";
import { ReactionResultCard } from "./reaction-result-card.jsx";

const MAX_TEXT_LENGTH = 200;
const TEXT_FILTERS = [
  "workgroup",
  "workbook",
  "keyword",
  "reagent",
  "dateFrom",
  "dateTo",
];
// at least one of these is needed; workgroup and workbook only narrow a search
const SEARCH_TERMS = ["keyword", "reagent", "dateFrom", "dateTo", "structure"];

const INITIAL_VALUES = {
  workgroup: "",
  workbook: "",
  keyword: "",
  reagent: "",
  dateFrom: "",
  dateTo: "",
  structureMode: "substructure",
  roles: [...STRUCTURE_SEARCH_ROLES],
};

const validationSchema = object({
  keyword: string().max(
    MAX_TEXT_LENGTH,
    `Keep the keyword under ${MAX_TEXT_LENGTH} characters`,
  ),
  reagent: string().max(
    MAX_TEXT_LENGTH,
    `Keep the name under ${MAX_TEXT_LENGTH} characters`,
  ),
  dateFrom: string(),
  // YYYY-MM-DD strings compare correctly as text
  dateTo: string().test(
    "after-start",
    "The end date must be on or after the start date",
    (value, { parent }) =>
      !value || !parent.dateFrom || value >= parent.dateFrom,
  ),
  roles: array().of(string()),
});

/**
 * Builds the search API request from the form and the drawn structure. Empty
 * fields are left out.
 * @param {object} values The form values.
 * @param {string} smiles The drawn structure, or "".
 * @returns {object} The request body.
 */
const buildCriteria = (values, smiles) => {
  const criteria = {};
  for (const key of TEXT_FILTERS) {
    const value = values[key].trim();
    if (value) criteria[key] = value;
  }
  if (smiles) {
    criteria.structure = {
      smiles,
      mode: values.structureMode,
      roles: values.roles,
    };
  }
  return criteria;
};

/** Heading at the top of a panel. */
const PanelHeading = (props) => (
  <Heading as="h2" m={0} fontSize="18px" fontWeight={600} {...props} />
);

/**
 * Workgroup, workbook, keyword, reagent and date filters.
 * @param {object} props
 * @param {object} props.values The form values.
 * @param {(field: string, value: any) => void} props.setFieldValue Formik setter.
 * @param {{ name: string, workbooks: string[] }[]} props.workgroups The user's
 *   workgroups.
 */
const FilterPanel = ({ values, setFieldValue, workgroups }) => {
  const workbooks =
    workgroups.find((wg) => wg.name === values.workgroup)?.workbooks ?? [];

  return (
    <Card as="section" aria-labelledby="search-filters-heading" p="20px">
      <PanelHeading id="search-filters-heading" mb="16px">
        Filters
      </PanelHeading>
      <VStack spacing="16px" align="stretch">
        <SelectField
          name="workgroup"
          label="Workgroup"
          onChange={(e) => {
            setFieldValue("workgroup", e.target.value);
            // workbooks belong to a workgroup, so choose again
            setFieldValue("workbook", "");
          }}
        >
          <option value="">All workgroups</option>
          {workgroups.map((wg) => (
            <option key={wg.name} value={wg.name}>
              {wg.name}
            </option>
          ))}
        </SelectField>
        <SelectField
          name="workbook"
          label="Workbook"
          isDisabled={!values.workgroup}
          helpText={values.workgroup ? undefined : "Choose a workgroup first."}
        >
          <option value="">All workbooks</option>
          {workbooks.map((wb) => (
            <option key={wb} value={wb}>
              {wb}
            </option>
          ))}
        </SelectField>
        <TextField
          name="keyword"
          label="Keyword"
          type="search"
          placeholder="Name, description or reaction ID"
        />
        <TextField
          name="reagent"
          label="Reagent or solvent name"
          type="search"
          placeholder="e.g. triethylamine"
        />
        <TextField name="dateFrom" label="Created from" type="date" />
        <TextField name="dateTo" label="Created to" type="date" />
      </VStack>
    </Card>
  );
};

/**
 * Sketcher, match type and the parts of a reaction to match the structure in.
 * @param {object} props
 * @param {object} props.values The form values.
 * @param {(field: string, value: any) => void} props.setFieldValue Formik setter.
 * @param {React.Ref} props.sketcherRef Ref for the sketcher.
 * @param {() => void} props.onReady Called when the sketcher has loaded.
 */
const StructurePanel = ({ values, setFieldValue, sketcherRef, onReady }) => {
  const isSubstructure = values.structureMode === "substructure";
  return (
    <Card as="section" aria-labelledby="search-structure-heading" p="20px">
      <PanelHeading id="search-structure-heading" mb="4px">
        Structure
      </PanelHeading>
      <Muted mb="16px" fontSize="14px">
        Optional. Leave the sketcher empty to search with the filters only.
      </Muted>

      <Flex wrap="wrap" gap="16px" align="flex-end" mb="16px">
        <Flex role="group" aria-label="Match type">
          <SegmentButton
            position="first"
            isPressed={isSubstructure}
            onClick={() => setFieldValue("structureMode", "substructure")}
          >
            {STRUCTURE_MODE_LABELS.substructure}
          </SegmentButton>
          <SegmentButton
            position="last"
            isPressed={!isSubstructure}
            borderLeftWidth={isSubstructure ? "1px" : 0}
            borderLeftColor={isSubstructure ? INV.accent : undefined}
            onClick={() => setFieldValue("structureMode", "exact")}
          >
            {STRUCTURE_MODE_LABELS.exact}
          </SegmentButton>
        </Flex>

        <FormControl as="fieldset" w="auto">
          <FormLabel
            as="legend"
            m={0}
            mb="4px"
            fontSize="12px"
            fontWeight={500}
            color={INV.muted}
          >
            Match in
          </FormLabel>
          <CheckboxGroup
            value={values.roles}
            onChange={(roles) => setFieldValue("roles", roles)}
          >
            <Stack
              direction="row"
              spacing="16px"
              wrap="wrap"
              minH="44px"
              align="center"
            >
              {STRUCTURE_SEARCH_ROLES.map((role) => (
                <Checkbox key={role} value={role}>
                  {ROLE_LABELS[role]}
                </Checkbox>
              ))}
            </Stack>
          </CheckboxGroup>
        </FormControl>
      </Flex>

      <Box
        border="1px solid"
        borderColor={INV.line}
        borderRadius="8px"
        overflow="hidden"
      >
        <StructureSketcher
          ref={sketcherRef}
          title="Structure to search for"
          onReady={onReady}
        />
      </Box>
    </Card>
  );
};

/**
 * The search results. Nothing is shown before the first search.
 * @param {object} props
 * @param {{ results: object[], truncated: boolean } | null} props.response
 *   The latest search response, or null before the first search.
 */
const Results = ({ response }) => {
  if (!response) return null;
  const { results, truncated } = response;
  const count = results.length;
  return (
    <Box as="section" aria-labelledby="search-results-heading">
      <PanelHeading id="search-results-heading" mb="12px" aria-live="polite">
        {count === 1 ? "1 reaction found" : `${count} reactions found`}
      </PanelHeading>
      {truncated && (
        <Muted mb="12px" fontSize="14px">
          Only the first {count} matches are shown. Add more filters to narrow
          the search.
        </Muted>
      )}
      {count === 0 ? (
        <Card>
          <EmptyState>No reactions match your search.</EmptyState>
        </Card>
      ) : (
        <VStack spacing="12px" align="stretch">
          {results.map((reaction) => (
            <ReactionResultCard key={reaction.id} reaction={reaction} />
          ))}
        </VStack>
      )}
    </Box>
  );
};

/**
 * Reaction search, routed at /search. Combines optional filters (workgroup,
 * workbook, keyword, reagent or solvent name, creation dates) with an optional
 * drawn structure, matched exactly or as a substructure in the chosen parts of
 * each reaction. Reuses the inventory pages' styling.
 */
export const SearchPage = () => {
  const { search } = useBackendApi();
  const { data } = useSearchWorkgroups();
  const workgroups = data.workgroups;
  const sketcherRef = useRef(null);
  const [isSketcherReady, setIsSketcherReady] = useState(false);
  const [response, setResponse] = useState(null);
  const [error, setError] = useState(null);

  const handleSubmit = async (values) => {
    setError(null);
    const smiles = isSketcherReady
      ? (await sketcherRef.current.getSmiles()).trim()
      : "";
    if (smiles && values.roles.length === 0) {
      setError("Choose at least one place to match the structure.");
      return;
    }
    const criteria = buildCriteria(values, smiles);
    if (!SEARCH_TERMS.some((key) => key in criteria)) {
      setError("Enter a keyword, reagent, date or structure to search for.");
      return;
    }
    try {
      setResponse(await search.searchReactions(criteria));
    } catch (e) {
      console.error(e);
      setError(await describeError(e));
    }
  };

  const handleClear = async (resetForm) => {
    resetForm();
    setError(null);
    setResponse(null);
    if (isSketcherReady) await sketcherRef.current.clear();
  };

  if (workgroups.length === 0) {
    return (
      <InvPage breadcrumb="Reactions" title="Search">
        <Card>
          <EmptyState>
            You must belong to a workgroup to search reactions.
          </EmptyState>
        </Card>
      </InvPage>
    );
  }

  return (
    <InvPage breadcrumb="Reactions" title="Search">
      <Formik
        initialValues={INITIAL_VALUES}
        validationSchema={validationSchema}
        onSubmit={handleSubmit}
      >
        {({ values, setFieldValue, resetForm, isSubmitting }) => (
          <Form noValidate>
            <Grid
              templateColumns={{ base: "1fr", lg: "300px minmax(0, 1fr)" }}
              gap="24px"
              alignItems="start"
            >
              <FilterPanel
                values={values}
                setFieldValue={setFieldValue}
                workgroups={workgroups}
              />
              <VStack spacing="16px" align="stretch" minW={0}>
                <StructurePanel
                  values={values}
                  setFieldValue={setFieldValue}
                  sketcherRef={sketcherRef}
                  onReady={() => setIsSketcherReady(true)}
                />
                {error && <ErrorMessage mb={0}>{error}</ErrorMessage>}
                <Flex gap="12px" justify="flex-end">
                  <InvButton
                    onClick={() => handleClear(resetForm)}
                    isDisabled={isSubmitting}
                  >
                    Clear
                  </InvButton>
                  <InvButton primary type="submit" isLoading={isSubmitting}>
                    Search
                  </InvButton>
                </Flex>
                <Results response={response} />
              </VStack>
            </Grid>
          </Form>
        )}
      </Formik>
    </InvPage>
  );
};
