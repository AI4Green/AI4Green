import { Flex, Text, VStack } from "@chakra-ui/react";
import { StructureSketcher } from "components/core/structure-sketcher/structure-sketcher.jsx";
import { STRUCTURE_MODE_LABELS } from "components/core/structure-sketcher/structure-search-labels.js";
import { useBackendApi } from "contexts";
import { useRef, useState } from "react";
import { SegmentButton } from "./inventory-filters.jsx";
import { ErrorMessage, INV, InvModal } from "./ui.jsx";
import { describeError } from "./use-inventory-action.js";

/**
 * Draw a structure and find the workgroup's compounds that match it, exactly
 * or as a substructure. Only reads data, so any member can use it.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup whose inventory is searched.
 * @param {"exact" | "substructure"} [props.initialMode] The mode to start with.
 * @param {(match: { ids: Set<number>, mode: string }) => void} props.onApply
 *   Called with the matching compound ids.
 * @param {() => void} props.onClose Closes the modal.
 */
export const StructureSearchModal = ({
  workgroupName,
  initialMode = "substructure",
  onApply,
  onClose,
}) => {
  const { inventory } = useBackendApi();
  const sketcherRef = useRef(null);
  const [isReady, setIsReady] = useState(false);
  const [mode, setMode] = useState(initialMode);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);

  const handleSearch = async () => {
    setError(null);
    const smiles = await sketcherRef.current.getSmiles();
    if (!smiles.trim()) {
      setError("Draw a structure to search for.");
      return;
    }
    setIsLoading(true);
    try {
      const { workgroupCompoundIds } = await inventory.structureSearch(
        workgroupName,
        { smiles, mode },
      );
      onApply({ ids: new Set(workgroupCompoundIds), mode });
      onClose();
    } catch (e) {
      console.error(e);
      setError(await describeError(e));
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <InvModal
      isOpen
      onClose={onClose}
      title="Search by structure"
      actionLabel="Search"
      onAction={isReady ? handleSearch : undefined}
      isLoading={isLoading || !isReady}
      size="4xl"
    >
      <VStack align="stretch" spacing="16px">
        {error && <ErrorMessage>{error}</ErrorMessage>}
        <Text color={INV.muted}>
          Draw a structure, then find the compounds in this inventory that are
          that structure or contain it.
        </Text>
        <Flex role="group" aria-label="Match type">
          <SegmentButton
            position="first"
            isPressed={mode === "substructure"}
            onClick={() => setMode("substructure")}
          >
            {STRUCTURE_MODE_LABELS.substructure}
          </SegmentButton>
          <SegmentButton
            position="last"
            isPressed={mode === "exact"}
            borderLeftWidth={mode === "substructure" ? "1px" : 0}
            borderLeftColor={mode === "substructure" ? INV.accent : undefined}
            onClick={() => setMode("exact")}
          >
            {STRUCTURE_MODE_LABELS.exact}
          </SegmentButton>
        </Flex>
        <StructureSketcher
          ref={sketcherRef}
          title="Structure to search for"
          onReady={() => setIsReady(true)}
        />
      </VStack>
    </InvModal>
  );
};
