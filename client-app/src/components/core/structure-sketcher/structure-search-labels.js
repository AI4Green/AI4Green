/*
 * Display labels for structure searches, shared by the reaction search page
 * and the inventory. Kept out of the component files so they only export
 * components (needed for React fast refresh).
 */

/** Labels for the structure search modes. */
export const STRUCTURE_MODE_LABELS = {
  exact: "Exact match",
  substructure: "Contains substructure",
};

/** Plural labels for the parts of a reaction, e.g. for checkboxes. */
export const ROLE_LABELS = {
  reactant: "Reactants",
  reagent: "Reagents",
  solvent: "Solvents",
  product: "Products",
};

/** Singular labels for the parts of a reaction, e.g. "Matched in solvent". */
export const SINGULAR_ROLE_LABELS = {
  reactant: "reactant",
  reagent: "reagent",
  solvent: "solvent",
  product: "product",
};
