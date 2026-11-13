import { Flex } from "@chakra-ui/react";
import { useChemInventoryConnection, useWorkgroupInventory } from "api";
import { useMemo, useState } from "react";
import { Link as RouterLink, useParams } from "react-router-dom";
import InventoryBrowser from "./inventory-browser.jsx";
import { SegmentButton } from "./inventory-filters.jsx";
import ChemInventoryView from "./cheminventory-view.jsx";
import { ChemInventoryModal } from "./cheminventory-modal.jsx";
import { AddStockModal } from "./add-stock-modal.jsx";
import { HistoryDrawer } from "./history-drawer.jsx";
import { ImportModal } from "./import-modal.jsx";
import { AdjustModal, DisposeModal, MoveModal } from "./item-modals.jsx";
import { INV, InvButton, InvPage } from "./ui.jsx";

const SOURCES = {
  ai4green: "ai4green",
  chemInventory: "cheminventory",
};

/**
 * Switches between the AI4Green inventory and the workgroup's ChemInventory.
 * @param {object} props
 * @param {string} props.source The selected source.
 * @param {(source: string) => void} props.onChange Called with the new source.
 */
const SourceToggle = ({ source, onChange }) => (
  <Flex role="group" aria-label="Inventory source">
    <SegmentButton
      position="first"
      isPressed={source === SOURCES.ai4green}
      onClick={() => onChange(SOURCES.ai4green)}
    >
      AI4Green
    </SegmentButton>
    <SegmentButton
      position="last"
      isPressed={source === SOURCES.chemInventory}
      borderLeftWidth={source === SOURCES.ai4green ? "1px" : 0}
      borderLeftColor={source === SOURCES.ai4green ? INV.accent : undefined}
      onClick={() => onChange(SOURCES.chemInventory)}
    >
      ChemInventory
    </SegmentButton>
  </Flex>
);

/**
 * Workgroup inventory page, routed at /inventory/:workgroupName.
 * Shows the AI4Green inventory, or, when the server has ChemInventory set up,
 * a search of the workgroup's ChemInventory stock.
 * Loading the AI4Green inventory is handled by the surrounding Suspense
 * boundary; the ChemInventory connection loads separately so it can't hold
 * the page up.
 */
export default function InventoryPage() {
  const { workgroupName } = useParams();
  const { data: inventory } = useWorkgroupInventory(workgroupName);
  const data = inventory.compounds;
  // Only principal investigators and senior researchers can change stock. The
  // backend enforces this too; here it just hides actions that would fail.
  const { canEdit } = inventory.permissions;

  const { data: connection, mutate: mutateConnection } =
    useChemInventoryConnection(workgroupName);
  const showSourceToggle = Boolean(connection?.available);
  const [source, setSource] = useState(SOURCES.ai4green);
  const isChemInventory = showSourceToggle && source === SOURCES.chemInventory;

  // The open modal or drawer, or null. `type` is an item action (adjust |
  // move | dispose | history), addContainer, addStock, import or chemInventory.
  const [dialog, setDialog] = useState(null);

  // The compound each container belongs to, so item actions can show its name
  const compoundByItemId = useMemo(
    () => new Map(data.flatMap((wc) => wc.items.map((item) => [item.id, wc]))),
    [data],
  );

  // action is one of: adjust | move | dispose | history
  const handleItemAction = (action, item) =>
    setDialog({ type: action, item, wc: compoundByItemId.get(item.id) });
  const handleAddContainer = (wc) => setDialog({ type: "addContainer", wc });
  const handleAddStock = () => setDialog({ type: "addStock" });
  const handleImport = () => setDialog({ type: "import" });
  const handleOpenChemInventorySettings = () =>
    setDialog({ type: "chemInventory" });
  const closeDialog = () => setDialog(null);

  const handleChemInventoryChanged = async (isConnected) => {
    await mutateConnection();
    if (isConnected) setSource(SOURCES.chemInventory);
  };

  // props shared by the modals and drawer for a single container
  const itemDialogProps = dialog?.item && {
    workgroupName,
    item: dialog.item,
    compound: dialog.wc.compound,
    onClose: closeDialog,
  };

  return (
    <InvPage
      breadcrumb={`Workgroups / ${workgroupName}`}
      title="Inventory"
      actions={
        isChemInventory ? (
          connection.canManage && (
            <InvButton onClick={handleOpenChemInventorySettings}>
              ChemInventory settings
            </InvButton>
          )
        ) : (
          <>
            <InvButton as={RouterLink} to="locations">
              Manage locations
            </InvButton>
            {canEdit && (
              <>
                <InvButton onClick={handleImport}>Import CSV</InvButton>
                <InvButton primary onClick={handleAddStock}>
                  + Add stock
                </InvButton>
              </>
            )}
          </>
        )
      }
    >
      {showSourceToggle && (
        <SourceToggle source={source} onChange={setSource} />
      )}

      {isChemInventory ? (
        <ChemInventoryView
          // a new connection may search a different inventory, so start afresh
          key={connection.connectedAt ?? "not-connected"}
          workgroupName={workgroupName}
          connection={connection}
          onOpenSettings={handleOpenChemInventorySettings}
        />
      ) : (
        <InventoryBrowser
          key="ai4green"
          compounds={data}
          canEdit={canEdit}
          onItemAction={handleItemAction}
          onAddContainer={handleAddContainer}
        />
      )}

      {dialog?.type === "adjust" && <AdjustModal {...itemDialogProps} />}
      {dialog?.type === "move" && (
        <MoveModal {...itemDialogProps} locations={inventory.locations} />
      )}
      {dialog?.type === "dispose" && <DisposeModal {...itemDialogProps} />}
      {dialog?.type === "history" && <HistoryDrawer {...itemDialogProps} />}
      {(dialog?.type === "addContainer" || dialog?.type === "addStock") && (
        <AddStockModal
          workgroupName={workgroupName}
          compounds={data}
          locations={inventory.locations}
          workbooks={inventory.workbooks}
          workgroupCompound={dialog.wc}
          onClose={closeDialog}
        />
      )}
      {dialog?.type === "import" && (
        <ImportModal
          workgroupName={workgroupName}
          workbooks={inventory.workbooks}
          onClose={closeDialog}
        />
      )}
      {dialog?.type === "chemInventory" && connection && (
        <ChemInventoryModal
          workgroupName={workgroupName}
          connection={connection}
          onChanged={handleChemInventoryChanged}
          onClose={closeDialog}
        />
      )}
    </InvPage>
  );
}
