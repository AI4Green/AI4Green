import { Button, HStack, Icon, Text } from "@chakra-ui/react";
// import { useIsProjectInstructor, useProject, useProjectGroupsList } from "api";
import { Breadcrumbs } from "components/core/breadcrumbs";
import { DataTable, DataTableGlobalFilter } from "components/core/data-table";
// import {
//   CreateOrEditModal,
//   DeleteModal,
//   LockProjectGroupNotesModal,
//   RemoveStudentModal,
//   StudentInviteModal,
// } from "components/project-group/modal";
// import { columns } from "components/project-group/table";
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
import { useNavigate, useParams, useSearchParams } from "react-router-dom";
import { buildProjectPath } from "routes/project";

export const ProjectOverview = () => {
  const { user } = useUser();
  const { projectId } = useParams();
  const navigate = useNavigate();

  const [searchValue, setSearchValue] = useState("");

  const { data: project } = useProject(projectId);
  const { data: isProjectInstructor } = useIsProjectInstructor(projectId);

  const { tableData, list, mutate } = useTableData(projectId);
  const [searchParams, setSearchParams] = useSearchParams();
  const action = searchParams.get("action");

  const breadcrumbItems = [
    { label: "Projects", href: "/projects" },
    {
      label: project?.name,
    },
  ];

  return (
    <DefaultContentLayout>
      <Breadcrumbs items={breadcrumbItems} />
      <HStack my={2} w="100%" justifyContent="space-between">
        <DefaultContentHeader
          header="Project Groups"
          icon={TITLE_ICON_COMPONENTS.ProjectGroup}
        />

        {user.permissions?.includes(
          PROJECTMANAGEMENT_PERMISSIONS.InviteInstructors,
        ) && (
          <Button
            onClick={() =>
              navigate(`/projects/${projectId}/instructors`, { replace: true })
            }
            leftIcon={
              <Icon
                as={TITLE_ICON_COMPONENTS.Instructors}
                fontSize={12}
                color="blue.500"
              />
            }
            size="xs"
            variant="outline"
            py={{ base: 3, md: 4 }}
          >
            <Text fontSize="xs" fontWeight="medium">
              Project Instructors
            </Text>
          </Button>
        )}
      </HStack>
      <DataTable data={tableData} globalFilter={searchValue} columns={columns}>
        <HStack flex={1} justifyContent="flex-start">
          <DataTableGlobalFilter
            searchValue={searchValue}
            setSearchValue={setSearchValue}
            placeholder="Search"
          />
          {user.permissions?.includes(
            PROJECTMANAGEMENT_PERMISSIONS.CreateProjectGroups,
          ) &&
            isProjectInstructor && (
              <NewButton onClick={() => setSearchParams({ action: "new" })} />
            )}
        </HStack>
      </DataTable>
      {isProjectInstructor && (
        <>
          {action === "new" && (
            <CreateOrEditModal list={list} mutate={mutate} />
          )}
          {action === "edit" && (
            <CreateOrEditModal list={list} mutate={mutate} />
          )}
          {action === "delete" && <DeleteModal list={list} mutate={mutate} />}
          {action === "remove-student" && (
            <RemoveStudentModal list={list} mutate={mutate} />
          )}
          {action === "invite-students" && (
            <StudentInviteModal list={list} mutate={mutate} />
          )}
          {action === "lock-notes" && (
            <LockProjectGroupNotesModal list={list} mutate={mutate} />
          )}
        </>
      )}
    </DefaultContentLayout>
  );
};

const useTableData = (id) => {
  const { data: list, mutate } = useProjectGroupsList(id);
  const tableData = useMemo(
    () =>
      list?.map((pg) => ({
        id: pg.id,
        name: pg.name,
        startDate: pg.startDate,
        planningDeadline: pg.planningDeadline,
        experimentDeadline: pg.experimentDeadline,
        project: pg.project,
        subRows: pg.students.map((student) => ({
          targetPath: buildProjectPath(id, pg.id, student.id),
          id: student.id,
          name: student.name,
          email: student.email,
        })),
      })),
    [list, id],
  );
  return { list, mutate, tableData: tableData ?? [] };
};
