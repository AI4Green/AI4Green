import { HStack, Text, Tooltip } from "@chakra-ui/react";
import { useProjectType } from "src/api/project-type.ts";
import { Badge } from "@/components/core/Badge.jsx";
import { Breadcrumbs } from "@/components/core/breadcrumbs.jsx";
import { Area } from "@/components/project-type/canvas/area.jsx";
import { Section } from "@/components/project-type/canvas/section";
import { TITLE_ICON_COMPONENTS } from "@/constants";
import { DefaultContentLayout } from "layouts/default";
import { useParams } from "react-router-dom";

export const ProjectTypeCanvas = () => {
  const { projectTypeId } = useParams();
  const { data: projectType } = useProjectType(projectTypeId);
  const breadcrumbs = [
    {
      label: "Home",
      href: "/",
    },
    {
      label: "Project Type Management",
      href: "/project-type-management",
    },
    {
      label: projectType.name,
    },
  ];
  return (
    <DefaultContentLayout>
      <Breadcrumbs items={breadcrumbs} />
      <HStack spacing={4}>
        <Tooltip
          label={projectType.description}
          hasArrow
          placement="right"
          fontSize="xs"
        >
          <Text fontWeight="medium">{projectType.name}</Text>
        </Tooltip>
        <Badge
          label="Project Type"
          colorScheme="gray"
          leftIcon={TITLE_ICON_COMPONENTS.ProjectType}
          fontSize="xxs"
        />
      </HStack>
      <Area />
      <HStack align="start" spacing={6} w="full">
        <Section projectType={projectType} />
      </HStack>
    </DefaultContentLayout>
  );
};
