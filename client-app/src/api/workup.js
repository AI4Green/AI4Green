import { useBackendApi } from "contexts";
import useSWR from "swr";

/*
 * Workups API (Webapp/sources/blueprints/api/workups), used by the workup and
 * purification cards on the reaction page.
 *
 * Workbook members can read; only the reaction's creator can change its steps,
 * and only until the reaction is locked (the backend returns 403 otherwise, and
 * `canEdit` in the steps response tells the UI). Every write sends a JSON body,
 * which a cross-site HTML form can't, so the backend refuses anything else.
 */

/** Builds the base path for one reaction, encoding each name for the URL. */
const reactionPath = ({ workgroup, workbook, reaction }) =>
  `workups/${encodeURIComponent(workgroup)}/${encodeURIComponent(
    workbook,
  )}/${encodeURIComponent(reaction)}`;

/** Builds the base path for one workbook. */
const workbookPath = ({ workgroup, workbook }) =>
  `workups/${encodeURIComponent(workgroup)}/${encodeURIComponent(workbook)}`;

/** Builds the path for a workbook's workup templates. */
const templatesPath = (workbookNames) =>
  `${workbookPath(workbookNames)}/templates`;

/**
 * Checks an id is a positive integer before it goes into a URL, so a bad value
 * can't point the request at a different API route.
 * @param {*} id A step or template id.
 * @returns {number} The id.
 * @throws {TypeError} If the id isn't a positive integer.
 */
const safeId = (id) => {
  if (!Number.isInteger(id) || id <= 0) {
    throw new TypeError(`Invalid workup id: ${id}`);
  }
  return id;
};

export const workupFetchKeys = {
  stepDefinitions: "workups/step-definitions",
  steps: (reactionNames, category) =>
    `${reactionPath(reactionNames)}?category=${encodeURIComponent(category)}`,
  templates: (workbookNames, category) =>
    `${templatesPath(workbookNames)}?category=${encodeURIComponent(category)}`,
  solvents: (workbookNames) => `${workbookPath(workbookNames)}/solvents`,
};

/**
 * Workup write actions, added to the backend API context as `workups`.
 * `reactionNames` is `{ workgroup, workbook, reaction }` (the reaction code).
 */
export const getWorkupApi = ({ api }) => ({
  /** Add an empty step to the end of a category. Resolves to the new step. */
  addStep: (reactionNames, category, stepType) =>
    api
      .post(`${reactionPath(reactionNames)}/steps`, {
        json: { category, stepType },
      })
      .json(),

  /**
   * Save a step, replacing its settings (`data`), notes, temperature and
   * materials. Resolves to the saved step.
   */
  updateStep: (
    reactionNames,
    stepId,
    { data, notes, temperature, materials },
  ) =>
    api
      .put(`${reactionPath(reactionNames)}/steps/${safeId(stepId)}`, {
        json: { data, notes, temperature, materials },
      })
      .json(),

  /**
   * Look a compound up by name or CAS number in the compound database and the
   * user's novel compounds. Resolves to `{ status: "found", compound }`,
   * `{ status: "not_found", suggestions }` or `{ status: "cas_not_found" }`.
   */
  lookupCompound: (workbookNames, query) =>
    api
      .get(`${workbookPath(workbookNames)}/compounds`, {
        searchParams: { q: query },
      })
      .json(),

  /** Put a category's steps in a new order (every step id, in order). */
  reorderSteps: (reactionNames, category, stepIds) =>
    api
      .put(`${reactionPath(reactionNames)}/steps/order`, {
        json: { category, stepIds: stepIds.map(safeId) },
      })
      .json(),

  deleteStep: (reactionNames, stepId) =>
    api.delete(`${reactionPath(reactionNames)}/steps/${safeId(stepId)}`),

  /** Append empty copies of a workbook template's steps. */
  applyTemplate: (reactionNames, templateId) =>
    api
      .post(`${reactionPath(reactionNames)}/apply-template`, {
        json: { templateId: safeId(templateId) },
      })
      .json(),

  /** Save a category's step types and order as a workbook template. */
  saveAsTemplate: (reactionNames, category, { name, description }) =>
    api
      .post(`${reactionPath(reactionNames)}/templates`, {
        json: { category, name, description },
      })
      .json(),

  /** Delete a template (only its creator can). */
  deleteTemplate: (workbookNames, templateId) =>
    api.delete(`${templatesPath(workbookNames)}/${safeId(templateId)}`),
});

/**
 * Get what each step type records: step type -> `{ fields, materials,
 * temperature }`. `fields` are settings (`{ key, label, kind, unit, choices,
 * maximum }`), `materials` are material slots (`{ slot, label, role, picker,
 * units, multiple, repeated, defaultFate }`) and `temperature` is the
 * temperature's label, or null if the step type has none.
 */
export const useStepDefinitions = () => {
  const { apiFetcher } = useBackendApi();
  return useSWR(workupFetchKeys.stepDefinitions, apiFetcher, {
    revalidateOnFocus: false,
  });
};

/**
 * Get a reaction's steps for one category: `{ category, canEdit, steps }`.
 * @param {{ workgroup: string, workbook: string, reaction: string }} reactionNames
 * @param {"WORKUP" | "PURIFICATION"} category
 */
export const useWorkupSteps = (reactionNames, category) => {
  const { apiFetcher } = useBackendApi();
  return useSWR(workupFetchKeys.steps(reactionNames, category), apiFetcher, {
    revalidateOnFocus: false,
  });
};

/**
 * Get a workbook's workup templates of one category. Each has `id`, `name`,
 * `description`, `steps` (step types in order), `creator` and `canDelete`.
 * @param {{ workgroup: string, workbook: string }} workbookNames
 * @param {"WORKUP" | "PURIFICATION"} category
 */
export const useWorkupTemplates = (workbookNames, category) => {
  const { apiFetcher } = useBackendApi();
  return useSWR(
    workupFetchKeys.templates(workbookNames, category),
    apiFetcher,
    {
      revalidateOnFocus: false,
    },
  );
};

/**
 * Get the workbook's solvent list for the coloured solvent picker:
 * `{ solvents: [{ name, flag, compoundId, novelCompoundId, solventName }],
 * colours }`, where `colours` are the user's CHEM21 colours.
 * @param {{ workgroup: string, workbook: string }} workbookNames
 */
export const useWorkupSolvents = (workbookNames) => {
  const { apiFetcher } = useBackendApi();
  return useSWR(workupFetchKeys.solvents(workbookNames), apiFetcher, {
    revalidateOnFocus: false,
  });
};
