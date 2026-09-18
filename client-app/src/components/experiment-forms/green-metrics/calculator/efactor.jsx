import { number, object } from "yup";

import { MetricsCalculator } from "./metrics";

/**
 * E-Factor Calculator
 */
export const EFactorCalculator = ({ name, isDisabled }) => {
  return (
    <MetricsCalculator
      name={name}
      title="E-Factor Calculator"
      fields={[
        { name: "wasteMass", label: "Total Mass of Waste Generated" },
        { name: "productMass", label: "Mass of Product Obtained" },
        { name: "eFactor", label: "E-Factor", isDisabled: true },
      ]}
      validationSchema={validationSchema}
      handleSubmit={handleSubmit}
      isDisabled={isDisabled}
    />
  );
};
const validationSchema = object().shape({
  wasteMass: number()
    .required("Waste Mass is required")
    .typeError("Waste Mass must be a numeric value"),
  productMass: number()
    .required("Product Mass is required")
    .typeError("Product Mass must be a numeric value")
    .moreThan(0, "Product Mass must be greater than 0"),
});

const handleSubmit = (values, helpers) => {
  const wasteMass = parseFloat(values.wasteMass);
  const productMass = parseFloat(values.productMass);

  const eFactor = wasteMass / productMass;

  helpers.setValue({ ...values, eFactor });
};
