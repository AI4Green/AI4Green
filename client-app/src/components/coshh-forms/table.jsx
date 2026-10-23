import { HStack } from "@chakra-ui/react";
import { useCoshhFormsList, useCoshhInstanceList } from "api/coshh-forms";
import { DataTable, DataTableGlobalFilter } from "components/core/data-table";
import { CoshhCreateModal } from "components/coshh-forms/form.jsx";
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
        <NewForm />
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

const NewForm = () => {
  const navigate = useNavigate();
  const [searchParams, setSearchParams] = useSearchParams();
  const isOpen = searchParams.get("action") === "new";

  const close = () => setSearchParams({});

  return (
    <>
      <NewButton onClick={() => setSearchParams({ action: "new" })} />
      <CoshhCreateModal
        isOpen={isOpen}
        onClose={close}
        onCreated={(data) => navigate(`${data.id}`)}
      />
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
        uuid: coshhInstance.uuid,
        reaction_code: coshhInstance.reaction_code,
        workbook: coshhInstance.workbook,
        workgroup: coshhInstance.workgroup,
        stage: coshhInstance.approval_status,
        template: coshhInstance.template_name,
        targetPath: `${coshhInstance.id}`,
      })),
    [coshhInstances],
  );

  return { data: tableData ?? [] };
};
