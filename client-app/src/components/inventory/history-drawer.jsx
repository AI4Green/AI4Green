import {
  Box,
  Drawer,
  DrawerBody,
  DrawerCloseButton,
  DrawerContent,
  DrawerHeader,
  DrawerOverlay,
  Flex,
  Spinner,
  Text,
} from "@chakra-ui/react";
import { useItemHistory } from "api";
import { useEffect, useState } from "react";
import { ItemSummary } from "./item-modals.jsx";
import {
  Badge,
  EmptyState,
  ErrorMessage,
  focusRing,
  INV,
  Muted,
} from "./ui.jsx";
import { describeError } from "./use-inventory-action.js";
import { formatQuantity } from "./utils.jsx";

/**
 * Labels and badge tones for the ledger's reasons, matching the backend's
 * StockChangeReason enum.
 */
const REASONS = {
  received: { label: "Received", tone: "ok" },
  used_in_reaction: { label: "Used in reaction", tone: "muted" },
  adjustment: { label: "Adjusted", tone: "low" },
  disposed: { label: "Disposed", tone: "exp" },
};

/** Formats a ledger timestamp as e.g. "03 Mar 2026, 14:05". */
const formatDateTime = (iso) =>
  new Date(iso).toLocaleString("en-GB", {
    day: "2-digit",
    month: "short",
    year: "numeric",
    hour: "2-digit",
    minute: "2-digit",
  });

/** Formats a signed ledger change, e.g. "+500 mL" or "−20 mL". */
const formatChange = (change, unit) => {
  const value = Number(change);
  // a real minus sign, so it reads correctly next to the plus
  const sign = value > 0 ? "+" : "−";
  return `${sign}${formatQuantity(Math.abs(value))} ${unit}`;
};

/**
 * One ledger entry: when, why, by whom, and how much the quantity changed.
 * @param {object} props
 * @param {object} props.entry The InventoryLedgerEntry.
 * @param {string} props.unit The container's unit, which every change is in.
 */
const HistoryEntry = ({ entry, unit }) => {
  const reason = REASONS[entry.reason] ?? {
    label: entry.reason,
    tone: "muted",
  };
  return (
    <Box as="li" py="12px" borderBottom="1px solid" borderColor={INV.rowLine}>
      <Flex justify="space-between" align="baseline" gap="12px">
        <Badge tone={reason.tone}>{reason.label}</Badge>
        <Text as="b" fontVariantNumeric="tabular-nums">
          {formatChange(entry.change, unit)}
        </Text>
      </Flex>
      <Muted mt="4px">
        {formatDateTime(entry.createdAt)}
        {entry.user?.name && ` by ${entry.user.name}`}
      </Muted>
      {entry.notes && (
        <Text mt="4px" fontSize="13px">
          {entry.notes}
        </Text>
      )}
    </Box>
  );
};

/**
 * The ledger list with its loading, error and empty states.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup the container belongs to.
 * @param {object} props.item The InventoryItem.
 */
const HistoryList = ({ workgroupName, item }) => {
  const { data, error, isLoading } = useItemHistory(workgroupName, item.id);
  const [errorMessage, setErrorMessage] = useState(null);

  useEffect(() => {
    if (!error) return;
    let cancelled = false;
    describeError(error).then((message) => {
      if (!cancelled) setErrorMessage(message);
    });
    return () => {
      cancelled = true;
    };
  }, [error]);

  if (error) {
    // describeError is async, so show a generic message until it resolves
    return (
      <ErrorMessage>
        {errorMessage ?? "Couldn't load the history."}
      </ErrorMessage>
    );
  }
  if (isLoading) {
    return (
      <Flex justify="center" py="32px">
        <Spinner color={INV.accent} label="Loading history" />
      </Flex>
    );
  }
  if (data.length === 0) {
    return <EmptyState>No changes recorded for this container.</EmptyState>;
  }

  // newest first
  const entries = [...data].reverse();
  return (
    <Box as="ol" listStyleType="none" m={0} p={0}>
      {entries.map((entry) => (
        <HistoryEntry key={entry.id} entry={entry} unit={item.unit} />
      ))}
    </Box>
  );
};

/**
 * Slide-out panel listing every quantity change to a container. Opens from
 * the right so the inventory table stays visible beside it.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup the container belongs to.
 * @param {object} props.item The InventoryItem.
 * @param {object} props.compound The item's compound.
 * @param {() => void} props.onClose Closes the drawer.
 */
export const HistoryDrawer = ({ workgroupName, item, compound, onClose }) => (
  <Drawer isOpen placement="right" size="sm" onClose={onClose}>
    <DrawerOverlay />
    <DrawerContent color={INV.ink}>
      <DrawerCloseButton
        w="44px"
        h="44px"
        top="8px"
        right="8px"
        _focusVisible={focusRing}
      />
      <DrawerHeader fontSize="18px" fontWeight={600} pr="56px">
        Container history
      </DrawerHeader>
      <DrawerBody fontSize="14px">
        <ItemSummary item={item} compound={compound} />
        <HistoryList workgroupName={workgroupName} item={item} />
      </DrawerBody>
    </DrawerContent>
  </Drawer>
);
