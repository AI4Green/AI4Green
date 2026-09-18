import { HStack, Text, Tooltip } from "@chakra-ui/react";
import { useCoshhForm } from "api/coshh-forms";
import { Badge } from "components/core/Badge";
import { Breadcrumbs } from "components/core/breadcrumbs";
import { Area } from "components/coshh-forms/canvas/area";
import { Section } from "components/coshh-forms/canvas/section";
import { TITLE_ICON_COMPONENTS } from "constants";
import { DefaultContentLayout } from "layouts/default";
import { useParams } from "react-router-dom";

export const WorkupCanvas = () => {
  // create mock workup until api is hooked up
  const workupMock = {
    name: "example workup template",
    description: "canvas config",
  };
  return (
    <DefaultContentLayout>
      <HStack spacing={4}>
        <Tooltip
          label={workupMock.description}
          hasArrow
          placement="right"
          fontSize="xs"
        >
          <Text fontWeight="medium">{workupMock.name}</Text>
        </Tooltip>
        <Badge
          label="Workup Procedure"
          colorScheme="gray"
          leftIcon={TITLE_ICON_COMPONENTS.ProjectType}
          fontSize="xxs"
        />
      </HStack>
      {/* TODO ADD WORKUP SPECIFIC CANVAS*/}
      {/*<Area />*/}
      {/*<HStack align="start" spacing={6} w="full">*/}
      {/*  <Section coshhForm={coshhForm} />*/}
      {/*</HStack>*/}
    </DefaultContentLayout>
  );
};
