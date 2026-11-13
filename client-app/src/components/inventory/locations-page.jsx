import {
  Link,
  Table,
  Tbody,
  Text,
  Thead,
  Tr,
  VisuallyHidden,
} from "@chakra-ui/react";
import { useWorkgroupInventory } from "api";
import { useMemo, useState } from "react";
import { Link as RouterLink, useParams } from "react-router-dom";
import { ArchiveLocationModal, LocationFormModal } from "./location-modals.jsx";
import {
  Badge,
  Card,
  EmptyState,
  INV,
  InvButton,
  InvPage,
  InvTd,
  InvTh,
  LinkButton,
  Muted,
  tableProps,
  underline,
} from "./ui.jsx";

/**
 * Manage a workgroup's storage locations, routed at
 * /inventory/:workgroupName/locations. Lists every location with how many
 * containers it holds, and lets members add, edit, archive and restore them.
 * Loading is handled by the surrounding Suspense boundary.
 */
export default function LocationsPage() {
  const { workgroupName } = useParams();
  const { data: inventory } = useWorkgroupInventory(workgroupName);
  const { locations, compounds } = inventory;
  // hides changes the backend would refuse (it enforces this too)
  const { canEdit } = inventory.permissions;

  // The open modal, or null. `type` is add, edit or archive.
  const [dialog, setDialog] = useState(null);
  const closeDialog = () => setDialog(null);

  const inventoryPath = `/inventory/${encodeURIComponent(workgroupName)}`;

  // Active (not archived) containers in each location, keyed by location id
  const containerCounts = useMemo(() => {
    const counts = new Map();
    compounds
      .flatMap((wc) => wc.items)
      .filter((item) => !item.isArchived)
      .forEach((item) => {
        const id = item.location.id;
        counts.set(id, (counts.get(id) ?? 0) + 1);
      });
    return counts;
  }, [compounds]);

  // active locations first, then alphabetical
  const sortedLocations = useMemo(
    () =>
      [...locations].sort(
        (a, b) =>
          Number(a.isArchived) - Number(b.isArchived) ||
          a.name.localeCompare(b.name),
      ),
    [locations],
  );

  return (
    <InvPage
      breadcrumb={
        <>
          Workgroups / {workgroupName} /{" "}
          <Link
            as={RouterLink}
            to={inventoryPath}
            color={INV.accent}
            _hover={underline}
            _focusVisible={underline}
          >
            Inventory
          </Link>
        </>
      }
      title="Locations"
      actions={
        <>
          <InvButton as={RouterLink} to={inventoryPath}>
            Back to inventory
          </InvButton>
          {canEdit && (
            <InvButton primary onClick={() => setDialog({ type: "add" })}>
              + Add location
            </InvButton>
          )}
        </>
      }
    >
      <Card overflowX="auto">
        {sortedLocations.length === 0 ? (
          <EmptyState>
            No locations yet. Add one to start recording where stock is kept.
          </EmptyState>
        ) : (
          <Table {...tableProps}>
            <Thead>
              <Tr>
                <InvTh>Location</InvTh>
                <InvTh>Containers</InvTh>
                <InvTh>Status</InvTh>
                <InvTh>
                  <VisuallyHidden>Actions</VisuallyHidden>
                </InvTh>
              </Tr>
            </Thead>
            <Tbody>
              {sortedLocations.map((location) => (
                <Tr key={location.id}>
                  <InvTd>
                    <Text as="b">{location.name}</Text>
                    {location.description && (
                      <Muted>{location.description}</Muted>
                    )}
                  </InvTd>
                  <InvTd>{containerCounts.get(location.id) ?? 0}</InvTd>
                  <InvTd>
                    <Badge tone={location.isArchived ? "muted" : "ok"}>
                      {location.isArchived ? "Archived" : "Active"}
                    </Badge>
                  </InvTd>
                  <InvTd whiteSpace="nowrap" textAlign="right">
                    {canEdit && (
                      <>
                        <LinkButton
                          aria-label={`Edit ${location.name}`}
                          onClick={() => setDialog({ type: "edit", location })}
                        >
                          Edit
                        </LinkButton>
                        <LinkButton
                          aria-label={`${
                            location.isArchived ? "Restore" : "Archive"
                          } ${location.name}`}
                          onClick={() =>
                            setDialog({ type: "archive", location })
                          }
                        >
                          {location.isArchived ? "Restore" : "Archive"}
                        </LinkButton>
                      </>
                    )}
                  </InvTd>
                </Tr>
              ))}
            </Tbody>
          </Table>
        )}
      </Card>

      {(dialog?.type === "add" || dialog?.type === "edit") && (
        <LocationFormModal
          workgroupName={workgroupName}
          locations={locations}
          location={dialog.location}
          onClose={closeDialog}
        />
      )}
      {dialog?.type === "archive" && (
        <ArchiveLocationModal
          workgroupName={workgroupName}
          location={dialog.location}
          containerCount={containerCounts.get(dialog.location.id) ?? 0}
          onClose={closeDialog}
        />
      )}
    </InvPage>
  );
}
