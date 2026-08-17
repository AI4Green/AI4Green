import { Flex, Icon, Text, useDisclosure, Spinner } from "@chakra-ui/react";
import { useProjectTypesList } from "api";
import { ActionButton } from "components/core/action-button";
import { DataTableColumnHeader } from "components/core/data-table";
import { DeleteModal } from "components/project-type/modal-delete";
import { CreateOrEditProjectTypeModal } from "components/project-type/modal-form";
import { MoveStageModal } from "components/stage/move-stage";
import {
  PROJECT_TYPE_MANAGEMENT_PERMISSIONS,
  STAGE_TYPES,
  STAGES,
  STAGES_PERMISSIONS,
  STATUS_ICON_COMPONENTS,
  TITLE_ICON_COMPONENTS,
} from "constants";
import { useUser } from "contexts";
import {
  FaCheck,
  FaDraftingCompass,
  FaLink,
  FaPaperPlane,
  FaTimes,
  FaTrash,
} from "react-icons/fa";
import { useSearchParams } from "react-router-dom";

export const columns = [
  {
    accessorKey: "name",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Name" />
    ),
    cell: ({ cell }) => (
      <Flex alignItems="center" gap={2}>
        <Icon
          as={TITLE_ICON_COMPONENTS.ProjectType}
          color="green.600"
          fontSize="lg"
        />
        <Text>{cell.getValue()}</Text>
      </Flex>
    ),
    meta: {
      width: 1,
    },
  },
  {
    accessorKey: "description",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Description" />
    ),
    meta: {
      width: "sm",
    },
  },
  {
    accessorKey: "stage",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Status" />
    ),
    cell: ({ row, cell }) => (
      <Flex alignItems="center" gap={2}>
        <Icon
          as={STATUS_ICON_COMPONENTS[row.original.stage].icon}
          color={STATUS_ICON_COMPONENTS[row.original.stage].color}
        />
        {cell.getValue()}
      </Flex>
    ),
    meta: {
      width: 1,
    },
  },
  {
    accessorKey: "inUseCount",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="In Use Count" />
    ),
    meta: {
      width: 1,
    },
  },
  {
    id: "inUse",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="In Use" />
    ),
    cell: ({ row }) => {
      const count = row.original.inUseCount;
      return (
        <Flex alignItems="center" gap={2}>
          <Icon
            as={count > 0 ? FaCheck : FaTimes}
            color={count > 0 ? "green.600" : "red.600"}
            fontSize="lg"
          />
          {count > 0 ? "Yes" : "No"}
        </Flex>
      );
    },

    meta: {
      width: 1,
    },
  },
  {
    id: "actions",
    header: ({ column }) => (
      <DataTableColumnHeader column={column} title="Actions" />
    ),
    cell: ({ row }) => (
      <Action projectType={row.original} inUse={row.original.inUseCount > 0} />
    ),
    meta: {
      width: 1,
    },
  },
];

const Action = ({ projectType, inUse }) => {
  const { user, isLoading } = useUser();
  if (isLoading) return <Spinner boxSize={16} />;

  const [searchParams, setSearchParams] = useSearchParams();
  const action = searchParams.get("action");
  const id = searchParams.get("id");

  const {
    isOpen: isPublishOpen,
    onOpen: onPublishOpen,
    onClose: onPublishClose,
  } = useDisclosure();

  const {
    isOpen: isDraftOpen,
    onOpen: onDraftOpen,
    onClose: onDraftClose,
  } = useDisclosure();

  const {
    isOpen: isDeprecateOpen,
    onOpen: onDeprecateOpen,
    onClose: onDeprecateClose,
  } = useDisclosure();

  const { mutate } = useProjectTypesList();
  const modalDefaultProps = {
    record: {
      id: projectType.id,
      title: projectType.name,
    },
    type: STAGE_TYPES.ProjectType,
    mutate: mutate,
  };

  const actions = {
    edit: {
      isEligible: () =>
        user.permissions.includes(
          PROJECT_TYPE_MANAGEMENT_PERMISSIONS.EditProjectTypes,
        ) && projectType.stage === STAGES.Draft,
      icon: <FaLink />,
      label: "Edit",
      onClick: () => setSearchParams({ action: "edit", id: projectType.id }),
    },
    delete: {
      isEligible: () =>
        !inUse &&
        user.permissions.includes(
          PROJECT_TYPE_MANAGEMENT_PERMISSIONS.DeleteProjectTypes,
        ) &&
        projectType.stage === STAGES.Draft,
      icon: <FaTrash />,
      label: "Delete",
      onClick: () => setSearchParams({ action: "delete", id: projectType.id }),
    },
    publish: {
      isEligible: () =>
        user.permissions.includes(
          PROJECT_TYPE_MANAGEMENT_PERMISSIONS.EditProjectTypes,
        ) && projectType.permissions.includes(STAGES_PERMISSIONS.CanPublish),

      icon: <FaPaperPlane />,
      label: "Publish",
      onClick: onPublishOpen,
    },
    draft: {
      isEligible: () =>
        user.permissions.includes(
          PROJECT_TYPE_MANAGEMENT_PERMISSIONS.EditProjectTypes,
        ) && projectType.permissions.includes(STAGES_PERMISSIONS.CanPutInDraft),

      icon: <FaDraftingCompass />,
      label: "Convert to Draft",
      onClick: onDraftOpen,
    },
    deprecate: {
      isEligible: () =>
        user.permissions.includes(
          PROJECT_TYPE_MANAGEMENT_PERMISSIONS.EditProjectTypes,
        ) && projectType.permissions.includes(STAGES_PERMISSIONS.CanDeprecate),
      icon: <FaTrash />,
      label: "Deprecate",
      onClick: onDeprecateOpen,
    },
  };
  return (
    <>
      <ActionButton actions={actions} size="xs" />
      {action === "edit" && id && <CreateOrEditProjectTypeModal />}
      {action === "delete" && id && <DeleteModal />}

      {isPublishOpen && (
        <MoveStageModal
          fixedNextStage={STAGES.Ready}
          modalTitle="Publish Project Type"
          modalMessage="Are you sure you want to publish this project type?"
          successMessage="Project type published successfully"
          isModalOpen={isPublishOpen}
          onModalClose={onPublishClose}
          {...modalDefaultProps}
        />
      )}

      {isDraftOpen && (
        <MoveStageModal
          fixedNextStage={STAGES.Draft}
          modalTitle="Convert to Draft"
          modalMessage="Are you sure you want to convert this project type to draft?"
          successMessage="Project type converted to draft successfully"
          isModalOpen={isDraftOpen}
          onModalClose={onDraftClose}
          {...modalDefaultProps}
        />
      )}

      {isDeprecateOpen && (
        <MoveStageModal
          fixedNextStage={STAGES.Deprecated}
          modalTitle="Deprecate Project Type"
          modalMessage="Are you sure you want to deprecate this project type?"
          successMessage="Project type deprecated successfully"
          isModalOpen={isDeprecateOpen}
          onModalClose={onDeprecateClose}
          {...modalDefaultProps}
        />
      )}
    </>
  );
};
