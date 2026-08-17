import { HStack } from "@chakra-ui/react";
// import { useProjectsList } from "api";
import { DataTable, DataTableGlobalFilter } from "components/core/data-table";
import { CreateOrEditModal, DeleteModal } from "components/project/modal";
import { columns } from "components/project/table";
import {
  PROJECTMANAGEMENT_PERMISSIONS,
  STAGES,
  TITLE_ICON_COMPONENTS,
} from "constants";
import { useUser } from "contexts";
import {
  DefaultContentHeader,
  DefaultContentLayout,
  NewButton,
} from "layouts/default";
import { useMemo, useState } from "react";
import { useSearchParams } from "react-router-dom";

export const ProjectList = () => {
  const { user } = useUser();
  const { tableData, list, mutate } = useTableData();
  const [searchValue, setSearchValue] = useState("");
  const [searchParams, setSearchParams] = useSearchParams();
  const action = searchParams.get("action");

  return (
    <DefaultContentLayout>
      <HStack>
        <DefaultContentHeader
          header="Projects"
          icon={TITLE_ICON_COMPONENTS.Project}
        />
      </HStack>
      <DataTable data={tableData} globalFilter={searchValue} columns={columns}>
        <HStack flex={1} spacing={4} justify="flex-start">
          <DataTableGlobalFilter
            searchValue={searchValue}
            setSearchValue={setSearchValue}
            placeholder="Search"
          />
          {user.permissions?.includes(
            PROJECTMANAGEMENT_PERMISSIONS.CreateProjects,
          ) && <NewButton onClick={() => setSearchParams({ action: "new" })} />}
        </HStack>
      </DataTable>
      {action === "new" && <CreateOrEditModal list={list} mutate={mutate} />}
      {action === "delete" && <DeleteModal list={list} mutate={mutate} />}
      {action === "edit" && <CreateOrEditModal list={list} mutate={mutate} />}
    </DefaultContentLayout>
  );
};

const useTableData = () => {
  const { data: list, mutate } = useProjectsList();
  const tableData = useMemo(
    () =>
      list?.map((project) => ({
        id: project.id,
        name: project.name,
        status: project.status || STAGES.OnGoing,
        projectType: project.projectType,
        targetPath: `/projects/${project.id}`,
      })),
    [list],
  );
  return { list, mutate, tableData: tableData ?? [] };
};
