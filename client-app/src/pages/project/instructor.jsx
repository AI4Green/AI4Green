import { HStack } from "@chakra-ui/react";
// import { useProject, useProjectInstructors } from "api";
import { Breadcrumbs } from "components/core/breadcrumbs";
import { DataTable, DataTableGlobalFilter } from "components/core/data-table";
import {
  InstructorInviteModal,
  RemoveInstructorModal,
} from "components/project/modal";
import { instructorColumns } from "components/project/table";
import {
  PROJECTMANAGEMENT_PERMISSIONS,
  TITLE_ICON_COMPONENTS,
} from "constants";
import { useUser } from "contexts";
import {
  DefaultContentHeader,
  DefaultContentLayout,
  NewButton,
} from "layouts/default";
import { useMemo, useState } from "react";
import { useParams, useSearchParams } from "react-router-dom";

export const ProjectInstructorList = () => {
  const { user } = useUser();
  const { projectId } = useParams();
  const { data: project } = useProject(projectId);
  const { tableData, list, mutate } = useTableData(projectId);
  const [searchValue, setSearchValue] = useState("");
  const [searchParams, setSearchParams] = useSearchParams();
  const action = searchParams.get("action");

  const breadcrumbItems = [
    { label: "Projects", href: "/projects" },
    {
      label: project?.name,
      href: `/projects/${project?.id}`,
    },
    {
      label: "Project Instructors",
    },
  ];

  return (
    <DefaultContentLayout>
      <Breadcrumbs items={breadcrumbItems} />
      <HStack>
        <DefaultContentHeader
          header="Instructors"
          icon={TITLE_ICON_COMPONENTS.Instructors}
        />
      </HStack>
      <DataTable
        data={tableData}
        globalFilter={searchValue}
        columns={instructorColumns}
      >
        <HStack flex={1} spacing={4} justify="flex-start">
          <DataTableGlobalFilter
            searchValue={searchValue}
            setSearchValue={setSearchValue}
            placeholder="Search"
          />
          {user.permissions?.includes(
            PROJECTMANAGEMENT_PERMISSIONS.InviteInstructors,
          ) && (
            <NewButton
              onClick={() => setSearchParams({ action: "invite-instructors" })}
              label="Invite"
            />
          )}
        </HStack>
      </DataTable>
      {action === "invite-instructors" && (
        <InstructorInviteModal mutate={mutate} />
      )}
      {action === "remove-instructor" && (
        <RemoveInstructorModal list={list} mutate={mutate} />
      )}
    </DefaultContentLayout>
  );
};

const useTableData = (id) => {
  const { data: list, mutate } = useProjectInstructors(id);
  const tableData = useMemo(
    () =>
      list?.map((instructor) => ({
        id: instructor.id,
        name: instructor.fullName,
        email: instructor.email,
        roles: instructor.roles,
        emailConfirmed: instructor.emailConfirmed,
      })),
    [list],
  );
  return { list, mutate, tableData: tableData ?? [] };
};
