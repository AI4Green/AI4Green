import { HStack } from "@chakra-ui/react";
import { Breadcrumbs } from "components/core/breadcrumbs";
import {
  CoshhFormTable,
  CoshhInstanceTable,
} from "components/coshh-forms/table";
import { TITLE_ICON_COMPONENTS } from "constants/experiment-ui";
import {
  DefaultContentHeader,
  DefaultContentLayout,
} from "../../layouts/default";

export const CoshhFormsList = () => {
  const breadcrumbItems = [
    { label: "Home", href: "/" },
    {
      label: "COSHH Form Management",
    },
  ];
  return (
    <DefaultContentLayout>
      <Breadcrumbs items={breadcrumbItems} />
      <HStack>
        <DefaultContentHeader
          header="COSHH Form Management"
          icon={TITLE_ICON_COMPONENTS.ProjectType}
        />
      </HStack>
      <CoshhFormTable />
    </DefaultContentLayout>
  );
};

export const CoshhInstanceList = () => {
  const breadcrumbItems = [
    { label: "Home", href: "/" },
    {
      label: "My COSHH Forms",
    },
  ];
  return (
    <DefaultContentLayout>
      <Breadcrumbs items={breadcrumbItems} />
      <HStack>
        <DefaultContentHeader
          header=""
          icon={TITLE_ICON_COMPONENTS.ProjectType}
        />
      </HStack>
      <CoshhInstanceTable />
    </DefaultContentLayout>
  );
};
