import { HStack, Text, Tooltip } from "@chakra-ui/react";
import { useCoshhForm } from "api/project-type";
import { Badge } from "components/core/Badge";
import { Breadcrumbs } from "components/core/breadcrumbs";
import { Area } from "components/project-type/canvas/area";
import { Section } from "components/project-type/canvas/section";
import { TITLE_ICON_COMPONENTS } from "constants";
import { DefaultContentLayout } from "layouts/default";
import { useParams } from "react-router-dom";

export const CoshhFormCanvas = () => {
  const { coshhFormId } = useParams();
  const { data: coshhForm } = useCoshhForm(coshhFormId);
  const breadcrumbs = [
    {
      label: "Home",
      href: "/",
    },
    {
      label: "COSHH Form Management",
      href: "/coshh-form-management",
    },
    {
      label: coshhForm.name,
    },
  ];
  return (
    <DefaultContentLayout>
      <Breadcrumbs items={breadcrumbs} />
      <HStack spacing={4}>
        <Tooltip
          label={coshhForm.description}
          hasArrow
          placement="right"
          fontSize="xs"
        >
          <Text fontWeight="medium">{coshhForm.name}</Text>
        </Tooltip>
        <Badge
          label="COSHH Form"
          colorScheme="gray"
          leftIcon={TITLE_ICON_COMPONENTS.ProjectType}
          fontSize="xxs"
        />
      </HStack>
      <Area />
      <HStack align="start" spacing={6} w="full">
        <Section coshhForm={coshhForm} />
      </HStack>
    </DefaultContentLayout>
  );
};
