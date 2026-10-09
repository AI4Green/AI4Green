import { HStack } from "@chakra-ui/react";
import { useCoshhFormsList } from "api/coshh-forms";
import { DataTable, DataTableGlobalFilter } from "components/core/data-table";
import { columns } from "components/coshh-forms/columns";
import { CreateOrEditProjectTypeModal } from "components/coshh-forms/modal-form";
import { PROJECT_TYPE_MANAGEMENT_PERMISSIONS } from "constants";
import { useUser } from "contexts";
import { NewButton } from "layouts/default";
import { useMemo, useState } from "react";
import { useSearchParams } from "react-router-dom";

export const CoshhFormTable = () => {
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
  const { data: coshhForms } = useCoshhFormsList();
  const tableData = useMemo(
    () =>
      coshhForms?.map((coshhForm) => ({
        id: coshhForm.id,
        name: coshhForm.name,
        description: coshhForm.description,
        stage: coshhForm.stage,
        inUseCount: coshhForm.inUseCount,
        permissions: coshhForm.permissions,
        targetPath: `${coshhForm.id}`,
      })),
    [coshhForms],
  );

  return { data: tableData ?? [] };
};
