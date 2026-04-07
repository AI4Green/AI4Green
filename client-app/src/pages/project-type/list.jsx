import { HStack } from "@chakra-ui/react";
import { Breadcrumbs } from "components/core/breadcrumbs";
import { ProjectTypeTable } from "components/project-type/table";
import { TITLE_ICON_COMPONENTS } from "constants/experiment-ui";
import {
  DefaultContentHeader,
  DefaultContentLayout,
} from "../../layouts/default";

export const ProjectTypeList = () => {
  const breadcrumbItems = [
    { label: "Home", href: "/" },
    {
      label: "Project Type Management",
    },
  ];
  return (
    <DefaultContentLayout>
      <Breadcrumbs items={breadcrumbItems} />
      <HStack>
        <DefaultContentHeader
          header="Project Type Management"
          icon={TITLE_ICON_COMPONENTS.ProjectType}
        />
      </HStack>
      <ProjectTypeTable />
    </DefaultContentLayout>
  );
};
