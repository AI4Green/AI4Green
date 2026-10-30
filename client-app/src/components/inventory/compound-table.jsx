import {
  Box,
  Table,
  Tbody,
  Text,
  Thead,
  Tr,
  VisuallyHidden,
} from "@chakra-ui/react";
import { Fragment, useState } from "react";
import ContainersTable from "./containers-table.jsx";
import {
  Badge,
  BareIconButton,
  HazardTag,
  INV,
  InvTd,
  InvTh,
  LinkButton,
  Muted,
  tableProps,
} from "./ui.jsx";

const COLUMN_COUNT = 10;

/**
 * rows: [{ wc, summary, items }]
 *   wc       the WorkgroupCompound
 *   summary  result of summarise(wc)
 *   items    containers to show when the row is expanded (already filtered)
 * canEdit shows the actions that change stock; without it only History shows.
 */
export default function CompoundTable({
  rows,
  canEdit,
  onItemAction,
  onAddContainer,
}) {
  const [expanded, setExpanded] = useState(() => new Set());

  const toggle = (id) =>
    setExpanded((prev) => {
      const next = new Set(prev);
      next.has(id) ? next.delete(id) : next.add(id);
      return next;
    });

  return (
    <Table {...tableProps}>
      <Thead>
        <Tr>
          <InvTh w="44px">
            <VisuallyHidden>Expand</VisuallyHidden>
          </InvTh>
          <InvTh>Compound</InvTh>
          <InvTh>Hazards</InvTh>
          <InvTh>Total stock</InvTh>
          <InvTh>Containers</InvTh>
          <InvTh>Locations</InvTh>
          <InvTh>Status</InvTh>
          <InvTh>
            <VisuallyHidden>Actions</VisuallyHidden>
          </InvTh>
        </Tr>
      </Thead>
      <Tbody>
        {rows.map(({ wc, summary, items }) => {
          const isOpen = expanded.has(wc.id);
          const name = wc.compound.name;

          return (
            <Fragment key={wc.id}>
              <Tr bg={isOpen ? INV.rowOpen : undefined}>
                <InvTd>
                  <BareIconButton
                    aria-expanded={isOpen}
                    aria-label={`${isOpen ? "Collapse" : "Expand"} ${name}`}
                    onClick={() => toggle(wc.id)}
                    icon={
                      <svg
                        width="16"
                        height="16"
                        viewBox="0 0 16 16"
                        aria-hidden="true"
                      >
                        <polyline
                          points={isOpen ? "3,6 8,11 13,6" : "6,3 11,8 6,13"}
                          fill="none"
                          stroke="currentColor"
                          strokeWidth="2"
                        />
                      </svg>
                    }
                  />
                </InvTd>
                <InvTd>
                  <Text as="b">{name}</Text>
                  <Muted>{wc.compound.cas}</Muted>
                </InvTd>
                <InvTd>
                  {wc.compound.hazards.map((h) => (
                    <HazardTag key={h}>{h}</HazardTag>
                  ))}
                </InvTd>
                <InvTd>{summary.totalLabel}</InvTd>
                <InvTd>{summary.containerCount}</InvTd>
                <InvTd>{summary.locations.join(", ") || "-"}</InvTd>
                <InvTd>
                  <Badge tone={summary.status.tone}>
                    {summary.status.label}
                  </Badge>
                </InvTd>
                <InvTd>
                  {canEdit && (
                    <LinkButton
                      aria-label={`Add container of ${name}`}
                      onClick={() => onAddContainer(wc)}
                    >
                      Add container
                    </LinkButton>
                  )}
                </InvTd>
              </Tr>
              {isOpen && (
                <Tr>
                  <InvTd
                    colSpan={COLUMN_COUNT}
                    px={0}
                    py={0}
                    borderColor={INV.line}
                  >
                    <Box bg={INV.sub} pt={0} pr="14px" pb="12px" pl="56px">
                      <ContainersTable
                        items={items}
                        nested
                        canEdit={canEdit}
                        onItemAction={onItemAction}
                      />
                    </Box>
                  </InvTd>
                </Tr>
              )}
            </Fragment>
          );
        })}
      </Tbody>
    </Table>
  );
}
