import { number, object } from "yup";

import { MetricsCalculator } from "./metrics";

/**
 * Process Mass Intensity Calculator
 */
export const PMICalculator = ({ name, isDisabled }) => {
  return (
    <MetricsCalculator
      name={name}
      title="Process Mass Intensity Calculator"
      fields={[
        {
          name: "totalMassInProcess",
          label: "Total Mass in Process (including water):",
        },
        { name: "productMass", label: "Mass of Product Obtained" },
        { name: "pmi", label: "Process Mass Intensity", isDisabled: true },
      ]}
      validationSchema={validationSchema}
      handleSubmit={handleSubmit}
      isDisabled={isDisabled}
    />
  );
};
const validationSchema = object().shape({
  totalMassInProcess: number()
    .required("Total Mass in Process is required")
    .typeError("Total Mass in Process must be a numeric value"),
  productMass: number()
    .required("Product Mass is required")
    .typeError("Product Mass must be a numeric value")
    .moreThan(0, "Product Mass must be greater than 0"),
});

const handleSubmit = (values, helpers) => {
  const totalMassInProcess = parseFloat(values.totalMassInProcess);
  const productMass = parseFloat(values.productMass);

  const pmi = totalMassInProcess / productMass;

  helpers.setValue({ ...values, pmi });
};
