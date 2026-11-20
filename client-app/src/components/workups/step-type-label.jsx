import { HStack, Text } from "@chakra-ui/react";

import { getStepType } from "./step-types";

/**
 * A step type's icon and label, for menus and template previews.
 * @param {{ stepType: string }} props
 */
export const StepTypeLabel = ({ stepType }) => {
  const { icon: Icon, label } = getStepType(stepType);
  return (
    <HStack spacing={2}>
      <Icon aria-hidden />
      <Text as="span">{label}</Text>
    </HStack>
  );
};
