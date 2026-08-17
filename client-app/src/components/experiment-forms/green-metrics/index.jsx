import { FormLabel, VStack } from "@chakra-ui/react";
import { useField } from "formik";

import {
  EFactorCalculator,
  PMICalculator,
  RMECalculator,
  WasteIntensityCalculator,
} from "./calculator";

/**
 * Green Metrics Calculator.
 * It consists of Waste Intensity, E-Factor, and RME calculators
 */

export const GreenMetricsCalculator = ({ name, isDisabled }) => {
  useField(name);

  return (
    <VStack w="full" align="flex-start" p={6} spacing={8}>
      <FormLabel>Calculate Green Metrics</FormLabel>
      <WasteIntensityCalculator
        name={`${name}.wasteIntensityCalculation`}
        isDisabled={isDisabled}
      />

      <EFactorCalculator
        name={`${name}.efactorCalculation`}
        isDisabled={isDisabled}
      />

      <RMECalculator name={`${name}.rmeCalculation`} isDisabled={isDisabled} />

      <PMICalculator name={`${name}.pmiCalculation`} isDisabled={isDisabled} />
    </VStack>
  );
};
