export default function InventoryFilters({
  filters,
  onChange,
  locations,
  suppliers,
  view,
  onViewChange,
}) {
  const set = (key) => (e) =>
    onChange({
      ...filters,
      [key]: e.target.type === "checkbox" ? e.target.checked : e.target.value,
    });

  return (
    <div className="inv-filters">
      <div className="inv-field inv-field-grow">
        <label htmlFor="inv-search">Search</label>
        <input
          id="inv-search"
          type="search"
          placeholder="Name, CAS number, barcode or lot"
          value={filters.search}
          onChange={set("search")}
        />
      </div>

      <div className="inv-field">
        <label htmlFor="inv-location">Location</label>
        <select
          id="inv-location"
          value={filters.location}
          onChange={set("location")}
        >
          <option value="">All locations</option>
          {locations.map((l) => (
            <option key={l} value={l}>
              {l}
            </option>
          ))}
        </select>
      </div>

      <div className="inv-field">
        <label htmlFor="inv-supplier">Supplier</label>
        <select
          id="inv-supplier"
          value={filters.supplier}
          onChange={set("supplier")}
        >
          <option value="">All suppliers</option>
          {suppliers.map((s) => (
            <option key={s} value={s}>
              {s}
            </option>
          ))}
        </select>
      </div>

      <div className="inv-field">
        <label htmlFor="inv-status">Status</label>
        <select id="inv-status" value={filters.status} onChange={set("status")}>
          <option value="">Any status</option>
          <option value="low">Low stock</option>
          <option value="expiring">Expiring soon</option>
        </select>
      </div>

      <div className="inv-segmented" role="group" aria-label="View">
        <button
          type="button"
          aria-pressed={view === "compound"}
          onClick={() => onViewChange("compound")}
        >
          By compound
        </button>
        <button
          type="button"
          aria-pressed={view === "containers"}
          onClick={() => onViewChange("containers")}
        >
          All containers
        </button>
      </div>

      <label className="inv-check">
        <input
          type="checkbox"
          checked={filters.showArchived}
          onChange={set("showArchived")}
        />
        Show archived
      </label>
    </div>
  );
}
