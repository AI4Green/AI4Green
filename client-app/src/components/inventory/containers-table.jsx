import {
  Table,
  Tbody,
  Text,
  Thead,
  Tr,
  VisuallyHidden,
} from "@chakra-ui/react";
import {
  Badge,
  EmptyState,
  InvTd,
  InvTh,
  LinkButton,
  Muted,
  tableProps,
} from "./ui.jsx";
import { formatDate, formatQuantity, itemState } from "./utils.jsx";

const ACTIONS = [
  ["adjust", "Adjust"],
  ["move", "Move"],
  ["dispose", "Dispose"],
  ["history", "History"],
];

/**
 * One row per InventoryItem.
 * Pass showCompound for the flat "All containers" view, where items carry
 * compoundName and cas. Pass nested when rendered inside an expanded compound
 * row, which uses a denser layout.
 */
export default function ContainersTable({
  items,
  showCompound = false,
  nested = false,
  onItemAction,
}) {
  if (items.length === 0) {
    return <EmptyState>No containers to show.</EmptyState>;
  }

  return (
    <Table {...tableProps}>
      <Thead>
        <Tr>
          {showCompound && <InvTh nested={nested}>Compound</InvTh>}
          <InvTh nested={nested}>Location</InvTh>
          <InvTh nested={nested}>Quantity</InvTh>
          <InvTh nested={nested}>Lot</InvTh>
          <InvTh nested={nested}>Expiry</InvTh>
          <InvTh nested={nested}>State</InvTh>
          <InvTh nested={nested}>
            <VisuallyHidden>Actions</VisuallyHidden>
          </InvTh>
        </Tr>
      </Thead>
      <Tbody>
        {items.map((item) => {
          const state = item.isArchived
            ? { label: "Archived", tone: "muted" }
            : itemState(item);
          const actions = item.isArchived
            ? ACTIONS.filter(([key]) => key === "history")
            : ACTIONS;

          return (
            <Tr key={item.id}>
              {showCompound && (
                <InvTd nested={nested}>
                  <Text as="b">{item.compoundName}</Text>
                  <Muted>{item.cas}</Muted>
                </InvTd>
              )}
              <InvTd nested={nested}>{item.location.name}</InvTd>
              <InvTd nested={nested}>
                {formatQuantity(item.quantity)} {item.unit}
              </InvTd>
              <InvTd nested={nested}>{item.lotNumber || "-"}</InvTd>
              <InvTd nested={nested}>{formatDate(item.expiryDate)}</InvTd>
              <InvTd nested={nested}>
                <Badge tone={state.tone}>{state.label}</Badge>
              </InvTd>
              <InvTd nested={nested} whiteSpace="nowrap">
                {actions.map(([key, label]) => (
                  <LinkButton
                    key={key}
                    aria-label={`${label} container ${
                      item.lotNumber || item.id
                    }`}
                    onClick={() => onItemAction(key, item)}
                  >
                    {label}
                  </LinkButton>
                ))}
              </InvTd>
            </Tr>
          );
        })}
      </Tbody>
    </Table>
  );
}
