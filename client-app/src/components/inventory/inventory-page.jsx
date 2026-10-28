import { Box, Flex, Grid, Heading, HStack, Text } from "@chakra-ui/react";
import { useMemo, useState } from "react";
import { useParams } from "react-router-dom";
import { useWorkgroupInventory } from "api";
import { daysUntil, EXPIRY_WARNING_DAYS, summarise } from "./utils.jsx";
import InventoryFilters from "./inventory-filters.jsx";
import CompoundTable from "./compound-table.jsx";
import ContainersTable from "./containers-table.jsx";
import { Card, EmptyState, INV, InvButton, MOBILE, Muted } from "./ui.jsx";

/** Summary tile with a label and a large figure. */
const StatCard = ({ label, value, unit }) => (
  <Card px="20px" py="16px">
    <Muted>{label}</Muted>
    <Text fontSize="28px" fontWeight={600}>
      {value}{" "}
      {unit && (
        <Text as="span" fontSize="14px" fontWeight={400} color={INV.muted}>
          {unit}
        </Text>
      )}
    </Text>
  </Card>
);

const PAGE_SIZE = 10;
const INITIAL_FILTERS = {
  search: "",
  location: "",
  supplier: "",
  status: "",
  showArchived: false,
};

/**
 * Workgroup inventory page, routed at /inventory/:workgroupName.
 * Loading is handled by the surrounding Suspense boundary.
 */
export default function InventoryPage() {
  const { workgroupName } = useParams();
  const { data: inventory } = useWorkgroupInventory(workgroupName);
  const data = inventory.compounds;

  const [filters, setFilters] = useState(INITIAL_FILTERS);
  const [view, setView] = useState("compound"); // compound | containers
  const [page, setPage] = useState(0);

  const updateFilters = (next) => {
    setFilters(next);
    setPage(0);
  };
  const updateView = (next) => {
    setView(next);
    setPage(0);
  };

  // Options for the filter dropdowns
  const { locations, suppliers } = useMemo(() => {
    const items = data.flatMap((wc) => wc.items);
    const uniq = (values) => [...new Set(values.filter(Boolean))].sort();
    return {
      locations: uniq(items.map((i) => i.location.name)),
      suppliers: uniq(items.map((i) => i.supplier)),
    };
  }, [data]);

  // Headline numbers, always for the whole workgroup (not the filtered view)
  const stats = useMemo(() => {
    const summaries = data.map(summarise);
    return {
      containers: summaries.reduce((n, s) => n + s.containerCount, 0),
      lowStock: summaries.filter((s) => s.belowMin).length,
      expiring: summaries.reduce((n, s) => n + s.expiring + s.expired, 0),
    };
  }, [data]);

  // Apply search and filters. Totals come from the full compound; only the
  // containers shown in the expanded row are narrowed by item-level filters.
  const rows = useMemo(() => {
    const q = filters.search.trim().toLowerCase();
    const itemFilterOn = Boolean(filters.location || filters.supplier);

    return data
      .map((wc) => ({
        wc,
        summary: summarise(wc),
        items: wc.items.filter(
          (i) =>
            (filters.showArchived || !i.isArchived) &&
            (!filters.location || i.location.name === filters.location) &&
            (!filters.supplier || i.supplier === filters.supplier),
        ),
      }))
      .filter(({ wc, summary, items }) => {
        if (itemFilterOn && items.length === 0) return false;
        if (filters.status === "low" && !summary.belowMin) return false;
        if (
          filters.status === "expiring" &&
          summary.expiring + summary.expired === 0
        ) {
          return false;
        }
        if (q) {
          const inCompound =
            wc.compound.name.toLowerCase().includes(q) ||
            wc.compound.cas.includes(q);
          const inItems = wc.items.some((i) =>
            [i.lotNumber, i.barcode].some(
              (v) => v && v.toLowerCase().includes(q),
            ),
          );
          if (!inCompound && !inItems) return false;
        }
        return true;
      });
  }, [data, filters]);

  const flatItems = useMemo(
    () =>
      rows.flatMap(({ wc, items }) =>
        items
          .filter(
            (i) => filters.status !== "expiring" || isExpiringOrExpired(i),
          )
          .map((i) => ({
            ...i,
            compoundName: wc.compound.name,
            cas: wc.compound.cas,
          })),
      ),
    [rows, filters.status],
  );

  const list = view === "compound" ? rows : flatItems;
  const total = list.length;
  const start = page * PAGE_SIZE;
  const pageList = list.slice(start, start + PAGE_SIZE);
  const noun = view === "compound" ? "compounds" : "containers";

  // Placeholders: wire these to modals and routes in the next pass.
  const handleItemAction = (action, item) => {
    // action is one of: adjust | move | dispose | history
    console.log("TODO open", action, "for item", item.id);
  };
  const handleAddContainer = (wc) =>
    console.log("TODO add container for", wc.id);
  const handleAddStock = () => console.log("TODO open add stock form");
  const handleManageLocations = () => console.log("TODO go to locations page");

  return (
    <Box bg={INV.bg} color={INV.ink} minH="100vh">
      <Flex
        direction="column"
        gap="24px"
        maxW="1280px"
        mx="auto"
        px="40px"
        py="32px"
        sx={{ [MOBILE]: { px: "16px", py: "20px" } }}
      >
        <Flex
          as="header"
          justify="space-between"
          align="flex-end"
          wrap="wrap"
          gap="16px"
        >
          <Box>
            <Text m={0} mb="4px" fontSize="12px" color={INV.muted}>
              Workgroups / {workgroupName}
            </Text>
            <Heading
              as="h1"
              m={0}
              fontSize="32px"
              fontWeight={600}
              lineHeight="inherit"
            >
              Inventory
            </Heading>
          </Box>
          <HStack spacing="12px">
            <InvButton onClick={handleManageLocations}>
              Manage locations
            </InvButton>
            <InvButton primary onClick={handleAddStock}>
              + Add stock
            </InvButton>
          </HStack>
        </Flex>

        <Grid
          as="section"
          aria-label="Summary"
          templateColumns="repeat(3, 1fr)"
          gap="16px"
          sx={{ [MOBILE]: { gridTemplateColumns: "1fr" } }}
        >
          <StatCard label="Containers in stock" value={stats.containers} />
          <StatCard label="Low stock" value={stats.lowStock} unit="compounds" />
          <StatCard
            label={`Expiring within ${EXPIRY_WARNING_DAYS} days`}
            value={stats.expiring}
            unit="containers"
          />
        </Grid>

        <InventoryFilters
          filters={filters}
          onChange={updateFilters}
          locations={locations}
          suppliers={suppliers}
          view={view}
          onViewChange={updateView}
        />

        <Card overflowX="auto">
          {total === 0 && (
            <EmptyState>No stock matches your search.</EmptyState>
          )}
          {total > 0 && (
            <>
              {view === "compound" ? (
                <CompoundTable
                  rows={pageList}
                  onItemAction={handleItemAction}
                  onAddContainer={handleAddContainer}
                />
              ) : (
                <ContainersTable
                  items={pageList}
                  showCompound
                  onItemAction={handleItemAction}
                />
              )}
              <Flex justify="space-between" align="center" px="14px" py="12px">
                <Muted as="span">
                  Showing {start + 1} to {Math.min(start + PAGE_SIZE, total)} of{" "}
                  {total} {noun}
                </Muted>
                <HStack spacing="8px">
                  <InvButton
                    isDisabled={page === 0}
                    onClick={() => setPage(page - 1)}
                  >
                    Previous
                  </InvButton>
                  <InvButton
                    isDisabled={start + PAGE_SIZE >= total}
                    onClick={() => setPage(page + 1)}
                  >
                    Next
                  </InvButton>
                </HStack>
              </Flex>
            </>
          )}
        </Card>
      </Flex>
    </Box>
  );
}

function isExpiringOrExpired(item) {
  const d = daysUntil(item.expiryDate);
  return d !== null && d <= EXPIRY_WARNING_DAYS;
}
