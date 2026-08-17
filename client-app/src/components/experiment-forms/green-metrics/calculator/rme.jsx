import { number, object } from "yup";

import { MetricsCalculator } from "./metrics";

/**
 * Reaction Mass Efficiency Calculator
 */
export const RMECalculator = ({ name, isDisabled }) => {
  return (
    <MetricsCalculator
      name={name}
      title="Reaction Mass Efficiency Calculator"
      fields={[
        { name: "productMass", label: "Mass of Product" },
        { name: "reactantMass", label: "Total Mass of Reactants used" },
        { name: "rme", label: "Reaction Mass Efficiency", isDisabled: true },
      ]}
      validationSchema={validationSchema}
      handleSubmit={handleSubmit}
      isDisabled={isDisabled}
    />
  );
};
const validationSchema = object().shape({
  productMass: number()
    .required("Mass of Product is required")
    .typeError("Mass of Product must be a numeric value"),
  reactantMass: number()
    .required("Total Mass of Reactants used is required")
    .typeError("Total Mass of Reactants used must be a numeric value")
    .moreThan(0, "Total Mass of Reactants used must be greater than 0"),
});

const handleSubmit = (values, helpers) => {
  const productMass = parseFloat(values.productMass);
  const reactantMass = parseFloat(values.reactantMass);

  const rme = (productMass / reactantMass) * 100;

  helpers.setValue({ ...values, rme });
};
