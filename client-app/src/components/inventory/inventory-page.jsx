import { useMemo, useState } from "react";
import { useParams } from "react-router-dom";
import { useWorkgroupInventory } from "api";
import { daysUntil, EXPIRY_WARNING_DAYS, summarise } from "./utils.jsx";
import InventoryFilters from "./inventory-filters.jsx";
import CompoundTable from "./compound-table.jsx";
import ContainersTable from "./containers-table.jsx";
// import "./inventory.css";

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
    <div className="inv-page">
      <div className="inv-inner">
        <header className="inv-header">
          <div>
            <p className="inv-crumb">Workgroups / {workgroupName}</p>
            <h1>Inventory</h1>
          </div>
          <div className="inv-header-actions">
            <button
              type="button"
              className="inv-btn"
              onClick={handleManageLocations}
            >
              Manage locations
            </button>
            <button
              type="button"
              className="inv-btn inv-btn-primary"
              onClick={handleAddStock}
            >
              + Add stock
            </button>
          </div>
        </header>

        <section className="inv-stats" aria-label="Summary">
          <div className="inv-card inv-stat">
            <div className="inv-cas">Containers in stock</div>
            <div className="inv-stat-value">{stats.containers}</div>
          </div>
          <div className="inv-card inv-stat">
            <div className="inv-cas">Low stock</div>
            <div className="inv-stat-value">
              {stats.lowStock} <span>compounds</span>
            </div>
          </div>
          <div className="inv-card inv-stat">
            <div className="inv-cas">
              Expiring within {EXPIRY_WARNING_DAYS} days
            </div>
            <div className="inv-stat-value">
              {stats.expiring} <span>containers</span>
            </div>
          </div>
        </section>

        <InventoryFilters
          filters={filters}
          onChange={updateFilters}
          locations={locations}
          suppliers={suppliers}
          view={view}
          onViewChange={updateView}
        />

        <div className="inv-card inv-table-card">
          {total === 0 && (
            <p className="inv-empty">No stock matches your search.</p>
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
              <div className="inv-pager">
                <span className="inv-cas">
                  Showing {start + 1} to {Math.min(start + PAGE_SIZE, total)} of{" "}
                  {total} {noun}
                </span>
                <div>
                  <button
                    type="button"
                    className="inv-btn"
                    disabled={page === 0}
                    onClick={() => setPage(page - 1)}
                  >
                    Previous
                  </button>
                  <button
                    type="button"
                    className="inv-btn"
                    disabled={start + PAGE_SIZE >= total}
                    onClick={() => setPage(page + 1)}
                  >
                    Next
                  </button>
                </div>
              </div>
            </>
          )}
        </div>
      </div>
    </div>
  );
}

function isExpiringOrExpired(item) {
  const d = daysUntil(item.expiryDate);
  return d !== null && d <= EXPIRY_WARNING_DAYS;
}
