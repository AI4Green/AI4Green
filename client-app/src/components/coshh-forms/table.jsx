import { HStack } from "@chakra-ui/react";
import { useCoshhFormsList, useCoshhInstanceList } from "api/coshh-forms";
import { DataTable, DataTableGlobalFilter } from "components/core/data-table";
import {
  coshhFormColumns,
  coshhInstanceColumns,
} from "components/coshh-forms/columns";
import { CreateOrEditProjectTypeModal } from "components/coshh-forms/modal-form";
import { PROJECT_TYPE_MANAGEMENT_PERMISSIONS } from "constants";
import { useUser } from "contexts";
import { NewButton } from "layouts/default";
import { useMemo, useState } from "react";
import { useSearchParams } from "react-router-dom";
import { useNavigate } from "react-router-dom";

export const CoshhFormTable = () => {
  const { user } = useUser();
  const { data } = useCoshhFormTableData();
  const navigate = useNavigate();
  const [searchValue, setSearchValue] = useState("");
  return (
    <DataTable
      data={data}
      columns={coshhFormColumns}
      globalFilter={searchValue}
      onRowClick={(row) => {
        if (row.targetPath) {
          navigate(row.targetPath);
        }
      }}
    >
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

export const CoshhInstanceTable = () => {
  const { user } = useUser();
  const { data } = useCoshhInstanceTableData();
  const navigate = useNavigate();
  const [searchValue, setSearchValue] = useState("");
  return (
    <DataTable
      data={data}
      columns={coshhInstanceColumns}
      globalFilter={searchValue}
      onRowClick={(row) => {
        if (row.targetPath) {
          navigate(row.targetPath);
        }
      }}
    >
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

const useCoshhFormTableData = () => {
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

const useCoshhInstanceTableData = () => {
  const { data: coshhInstances } = useCoshhInstanceList();
  const tableData = useMemo(
    () =>
      coshhInstances?.map((coshhInstance) => ({
        id: coshhInstance.id,
        name: coshhInstance.name,
        description: coshhInstance.description,
        stage: coshhInstance.stage,
        inUseCount: coshhInstance.inUseCount,
        permissions: coshhInstance.permissions,
        targetPath: `${coshhInstance.id}`,
      })),
    [coshhInstances],
  );

  return { data: tableData ?? [] };
};
