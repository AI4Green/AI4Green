import { number, object } from "yup";

import { MetricsCalculator } from "./metrics";

/**
 * Waste Intensity Calculator
 */
export const WasteIntensityCalculator = ({ name, isDisabled }) => {
  return (
    <MetricsCalculator
      name={name}
      title="Waste Intensity Calculator"
      fields={[
        { name: "waste", label: "Waste Produced" },
        { name: "output", label: "Productivity Output" },
        { name: "wasteIntensity", label: "Waste Intensity", isDisabled: true },
      ]}
      validationSchema={validationSchema}
      handleSubmit={handleSubmit}
      isDisabled={isDisabled}
    />
  );
};
const validationSchema = object().shape({
  waste: number()
    .required("Waste Produced is required")
    .typeError("Waste Produced must be a numeric value"),
  output: number()
    .required("Productivity Output is required")
    .typeError("Productivity Output must be a numeric value")
    .moreThan(0, "Productivity Output must be greater than 0"),
});

const handleSubmit = (values, helpers) => {
  const wasteValue = parseFloat(values.waste);
  const outputValue = parseFloat(values.output);

  const wasteIntensity = wasteValue / outputValue;

  helpers.setValue({ ...values, wasteIntensity });
};
