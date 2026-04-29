import { HStack } from "@chakra-ui/react";
import { useProjectTypesList } from "../../api/project-type.js";
import { DataTable, DataTableGlobalFilter } from "components/core/data-table";
import { columns } from "components/project-type/columns";
import { CreateOrEditProjectTypeModal } from "components/project-type/modal-form";
import { PROJECT_TYPE_MANAGEMENT_PERMISSIONS } from "constants";
import { useUser } from "contexts";
import { NewButton } from "layouts/default";
import { useMemo, useState } from "react";
import { useSearchParams } from "react-router-dom";

export const ProjectTypeTable = () => {
  // const { user } = useUser();
  const user = {
    id: "user-123",
    name: "Dev User",
    email: "dev@example.com",
    // Ensure this includes the specific permission the UI is looking for
    permissions: [
      "CreateProjectTypes",
      "EditProjectTypes",
      "DeleteProjectTypes",
    ],
  };
  const { data } = useTableData();
  const [searchValue, setSearchValue] = useState("");
  return (
    <DataTable data={data} columns={columns} globalFilter={searchValue}>
      <HStack flex={1} justifyContent="flex-start">
        <DataTableGlobalFilter
          searchValue={searchValue}
          setSearchValue={setSearchValue}
          placeholder="Search"
        />
        {user.permissions?.includes(
          PROJECT_TYPE_MANAGEMENT_PERMISSIONS.CreateProjectTypes,
        ) && <New />}
      </HStack>
    </DataTable>
  );
};

const New = () => {
  const [searchParams, setSearchParams] = useSearchParams();
  const action = searchParams.get("action");
  return (
    <>
      <NewButton onClick={() => setSearchParams({ action: "new" })} />
      {action === "new" && <CreateOrEditProjectTypeModal />}
    </>
  );
};

const useTableData = () => {
  const { data: projectTypes } = useProjectTypesList();
  const tableData = useMemo(
    () =>
      projectTypes?.map((projectType) => ({
        id: projectType.id,
        name: projectType.name,
        description: projectType.description,
        stage: projectType.stage,
        inUseCount: projectType.inUseCount,
        permissions: projectType.permissions,
        targetPath: `${projectType.id}`,
      })),
    [projectTypes],
  );

  return { data: tableData ?? [] };
};
